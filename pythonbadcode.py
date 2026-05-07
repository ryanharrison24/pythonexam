import sys

def validate_sequence(sequence, k):
    """
    Check if a DNA sequence is valid for kmer analysis.
    Returns False if the sequence is too short or contains numbers.
    Returns True if the sequence looks good.
    """
    # sequence has to be longer than k or we cant get a next character
    # fixed this from < to <= because a sequence exactly equal to k
    # has no character after the last kmer
    if len(sequence) <= k:
        return False
    
    # DNA should only have letters. Remove anything with numbers
    for nucleotide in sequence:
        if nucleotide in '1234567890':
            return False
    
    # passed both checks so its good
    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """
    Add a kmer and its next character to our running dictionary.
    If we havent seen this kmer before, create a new entry for it.
    If we have seen it, just increment the count.
    """
    # first time seeing this kmer so add it to the dict
    # fixed this to start at 0 so when we add 1 below it becomes 1
    # it was starting at 1 before so it was always one too high
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    # seen it before so add to count
    kmer_data[kmer]['count'] += 1
    
    # same for the next character
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    
    # increment the next character count
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k, kmer_data=None):
    """
    Go through a sequence and find all kmers of length k.
    For each kmer also record what character comes after it.
    Can take an existing kmer_data dict so multiple sequences
    can be combined together.
    Returns a dictionary with all the kmer data.
    """
    # if no existing data was passed in, start fresh
    # this lets us accumulate results across multiple sequences
    if kmer_data is None:
        kmer_data = {}
    
    # loop through the sequence stopping early enough that there is always a character after the kmer
    for i in range(len(sequence) - k):
        
        # pull out the kmer at this spot
        kmer = sequence[i:i+k]
        
        # and the character right after it
        next_char = sequence[i+k]
        
        # store it
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    return kmer_data

def write_results_to_file(kmer_data, output_filename):
    """
    Write all the kmer results out to a text file.
    Kmers are sorted alphabetically.
    Each line shows the kmer, its total count, and what characters follow it
    """
    # sort alphabetically so the output is easy to read
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            
            # get the count and next characters for this kmer
            count = kmer_data[kmer]['count']
            next_chars = kmer_data[kmer]['next_chars']
            
            # format it like "A:3 G:1 T:2"
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            # fixed this to also write the total count to the output
            f.write(f"{kmer} {count} {next_char_str}\n")

def main():
    """
    Main function that runs everything.
    Takes three arguments from the command line:
    - the input file with sequences
    - the value of k
    - the output file to write results to
    """
    # get the arguments from the command line
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")

    # fixed: kmer_data needs to be outside the loop so
    # all sequences get combined instead of overwriting each other
    kmer_data = {}

    # go through the input file line by line
    with open(sequence_file, 'r') as f:
        for sequence in f:
            
            # strip whitespace off the ends
            sequence = sequence.strip()

            # skip it if its not valid
            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            # fixed: pass kmer_data in so counts build up across all sequences
            kmer_data = count_kmers_with_context(sequence, k, kmer_data)
    
    # fixed: write the file once after all sequences are done
    write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
