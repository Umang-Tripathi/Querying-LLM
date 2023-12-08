from openai import OpenAI
import sys

MY_KEY = "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"
client = OpenAI(api_key=MY_KEY)


if __name__ == '__main__':
    
    # python3 <filename> API_KEY num_paragraphs query.txt
    if len(sys.argv) < 4:
        print("Usage: python3 api_call.py API_KEY num_paragraphs query.txt")
        sys.exit(1)

    # Read the API key from the command line
    
    num_paragraphs = int(sys.argv[2])
    

    # Read the paragraphs from the files
    paragraphs = []

    for i in range(num_paragraphs):
        filename = 'paragraph_' + str(i) + '.txt'
        
        with open(filename, 'r') as f:
            paragraphs.append(f.read())
            paragraphs.append('\n')
    
    # add query
    query_file = sys.argv[3]
    with open(query_file, 'r') as f:
        query = f.read()
        paragraphs.append(query)
        paragraphs.append('\n')

    # convert paragraphs to a single string
    paragraphs = '\n'.join(paragraphs)
    
    query = {
        "role": "user", "content": paragraphs
    }

    chat = client.chat.completions.create(model="gpt-3.5-turbo",
    messages=[query])

    reply = chat.choices[0].message.content

    print(reply)
