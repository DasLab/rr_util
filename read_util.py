import re

def get_total_mutdel(cols, ref_sequence):
    """
    Calculates the total number of mutations and deletions in a sequence alignment.

    Args:
        cols (list): A list of columns from a SAM file.
        ref_sequence (str): The reference sequence.

    Returns:
        tuple: (int: total mutations + deletions, str: aligned read sequence)
    """
    start_pos = int(cols[3])
    cigar = cols[5]
    signed_tmpl_len = int(cols[8])
    read = cols[9]

    seqa = ''
    read_pos = 0  # position in read
    seqa += '.' * (start_pos - 1)  # pad start

    # Parse CIGAR with regex: [('22', 'M'), ('1', 'I'), ...]
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    for length_str, op in cigar_tuples:
        length = int(length_str)
        if op == 'M':  # match or mismatch
            seqa += read[read_pos:read_pos + length]
            read_pos += length
        elif op == 'I':  # insertion to the reference — skip, not in alignment
            read_pos += length
        elif op == 'D':  # deletion from the reference — gap in read
            seqa += '-' * length
        elif op == 'S':  # soft clip — skip bases in read
            read_pos += length
        elif op == 'H':  # hard clip — ignore entirely
            continue
        elif op == 'N':  # skipped region — typically introns, same as D here
            seqa += '-' * length
        # P (padding) and others are rarely used; can be ignored or handled as needed

    # Ensure full alignment length
    if signed_tmpl_len < 0:
        seqa = '.' * (len(ref_sequence) - len(seqa)) + seqa
    else:
        seqa = seqa + '.' * (len(ref_sequence) - len(seqa))

    assert len(seqa) == len(ref_sequence)

    # Compare to reference for mutations and deletions
    pos = {'mut': [], 'del': []}
    for n, (ref_nt, read_nt) in enumerate(zip(ref_sequence, seqa)):
        if ref_nt not in 'ACGT-': continue
        if read_nt not in 'ACGT-': continue
        if read_nt == '-':
            pos['del'].append(n)
        elif ref_nt != read_nt:
            pos['mut'].append(n)

    return len(pos['mut']) + len(pos['del']), seqa

def get_md( cols ):
    for i,col in enumerate(cols):
        if len(col)>3 and col[:3] == 'MD:': return col[5:]
    return None

def parse_md(md):
    """ Parses the MD string and returns a list of tuples with matches and bases """
    tokens = re.findall(r'(\d+|\^[A-Za-z]+|[A-Za-z])', md)
    parsed_md = []
    for token in tokens:
        if token.isdigit():
            parsed_md.append((int(token), 'M'))  # Matches
        elif token.startswith('^'):
            parsed_md.append((len(token)-1, 'D', token[1:]))  # Deletions with sequence
        else:
            parsed_md.append((1, 'X', token))  # Single base mismatches
    return parsed_md

def get_align_read(start_pos, ref_seq, md):
    parsed_md = parse_md(md)
    ref_pos = start_pos - 1  # Convert 1-based start_pos to 0-based index
    align_read = list('-' * len(ref_seq))  # Initialize the alignment with gaps

    read_pos = 0
    for item in parsed_md:
        if item[1] == 'M':  # Matching bases
            for i in range(item[0]):
                if ref_pos < len(ref_seq):
                    align_read[ref_pos] = ref_seq[ref_pos]
                    ref_pos += 1
        elif item[1] == 'D':  # Deletion from reference
            ref_pos += item[0]
        elif item[1] == 'X':  # Mismatch
            if ref_pos < len(ref_seq):
                align_read[ref_pos] = item[2]
                ref_pos += 1

    return ''.join(align_read)

def parse_cigar(cigar):
    """
    Parse the CIGAR string into a list of tuples (operation, length).
    """
    return re.findall(r'(\d+)([MIDNSHP=X])', cigar)

def apply_cigar_to_read(reference, read, cigar, start_pos):
    """
    Apply the CIGAR operations to the reference and read to account for insertions and clipping.
    """
    ref_idx = start_pos - 1
    read_idx = 0
    read_with_cigar = []

    cigar_elements = parse_cigar(cigar)
    for length, op in cigar_elements:
        length = int(length)

        if op == 'M' or op == '=' or op == 'X':
            read_with_cigar.append(read[read_idx:read_idx + length])
            ref_idx += length
            read_idx += length
        elif op == 'I':
            #read_with_cigar.append(read[read_idx:read_idx + length])
            read_idx += length
        elif op == 'D':
            ref_idx += length
        elif op == 'N':
            ref_idx += length
        elif op == 'S':
            read_with_cigar.append(read[read_idx:read_idx + length])
            read_idx += length
        elif op == 'H':
            continue
        elif op == 'P':
            continue
    return ''.join(read_with_cigar)

def convert_md_string(reference, read, md, start_pos, cigar):
    """
    This function converts an MD string showing reference bases to one showing read bases.
    """
    read_from_cigar = apply_cigar_to_read(reference, read, cigar, start_pos)

    new_md = []
    ref_idx = start_pos - 1  # Convert 1-based start_pos to 0-based index
    read_idx = 0

    # Tokenize the MD string
    tokens = re.findall(r'(\d+|\^[A-Za-z]+|[A-Za-z])', md)

    for token in tokens:
        if token.isdigit():
            # Matches: copy the token directly
            ref_idx += int(token)
            read_idx += int(token)
            new_md.append(token)
        elif token.startswith('^'):
            # Deletions: copy token directly
            new_md.append(token)
            ref_idx += len(token) - 1  # Adjust ref_idx for the length of deleted bases
        else:
            # Mismatches: replace reference base with read base
            read_base = read_from_cigar[read_idx]
            new_md.append(read_base)
            ref_idx += 1
            read_idx += 1

    return ''.join(new_md)

def get_md_convert(cols,ref_seq):
    start_pos = int(cols[3])
    cigar = cols[5]
    read = cols[9]
    md = get_md(cols)
    md_convert = convert_md_string( ref_seq, read, md, start_pos, cigar )
    return str(start_pos),md_convert


def check_md_convert( start_md, ref_seq, align_read ):
    start_pos,md_convert = start_md
    infer_align_read = get_align_read( int(start_pos), ref_seq,  md_convert )
    if align_read.replace('.','-') != infer_align_read :
        print( start_pos, md_convert)
        print( ref_seq, 'ref' )
        print( align_read,'align' )
        print( infer_align_read, 'infer_align' )
        #print( read, 'read' )
        print()

