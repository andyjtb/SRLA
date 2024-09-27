import re
from deep_translator import GoogleTranslator
import os

def translate_comments(file_path):
    translator = GoogleTranslator(source='ja', target='en')
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Regular expression pattern to match comments
    pattern = r'//.*?$|/\*.*?\*/'
    matches = re.findall(pattern, content, re.MULTILINE | re.DOTALL)

    for match in matches:
        # Check if the comment contains Japanese characters
        if re.search(r'[^\x00-\xff]', match):
            try:
                translation = translator.translate(match)
                content = content.replace(match, f'{translation}' if '\n' not in match else f'/*\n{translation}\n*/')

                print(f'Translated comment in {file_path}:')
                print(f'Original: {match}')
                print(f'Translated: {translation}')
            except Exception as e:
                print(f'Error translating comment in {file_path}: {str(e)}')

    # Write the modified content back to the file
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)

def process_directory(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.c') or file.endswith('.h'):
                translate_comments(os.path.join(root, file))

process_directory(os.getcwd())