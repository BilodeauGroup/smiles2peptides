import re

def extract_characters(sequence):
    """
    Extracts individual characters from a sequence, considering anything inside
    curly braces `{}` as a single character.
    
    Args:
        sequence (str): The input sequence.
        
    Returns:
        list: A list of extracted characters, including groups inside curly braces as units.
    """
    # Regular expression to find groups inside `{}` or individual characters
    pattern = r'\{[^}]*\}|[a-zA-Z]'
    # Find all matches
    characters = re.findall(pattern, sequence)
    return characters
