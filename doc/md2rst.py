import requests

def md_to_rst(from_file, to_file):
    """
    @param from_file: {str} markdown
    @param to_file: {str} rst
    """
    response = requests.post(
        url='http://c.docverter.com/convert',
        data={'to': 'rst', 'from': 'markdown'},
        files={'input_files[]': open(from_file, 'rb')}
    )

    if response.ok:
        with open(to_file, "wb") as f:
            f.write(response.content)
    else:
        print("You need convert .md file to .rst online:")
        print("https://docverter.com/")

if __name__ == '__main__':
    md_to_rst("../README.md", "maindoc.rst")
