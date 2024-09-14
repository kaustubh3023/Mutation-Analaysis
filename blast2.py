import tkinter as tk
from tkinter import scrolledtext, messagebox
from Bio.Blast import NCBIWWW, NCBIXML
import requests
from bs4 import BeautifulSoup
import io
import traceback


def search_disorder(mutation):
    search_url = f"https://www.genecards.org/Search/Keyword?queryString={mutation}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'
    }
    
    try:
        print(f"Searching for mutation: {mutation}")
        response = requests.get(search_url, headers=headers)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        results = soup.find_all('div', class_='search-results')
        disorder_info = ""
        for result in results:
            title = result.find('a', class_='card-link')
            description = result.find('p', class_='description')
            if title and description:
                disorder_info += f"Title: {title.text.strip()}\nDescription: {description.text.strip()}\n\n"
                
        return disorder_info if disorder_info else "No disorder information found."

    except requests.exceptions.RequestException as e:
        return f"An error occurred while searching for disorders: {e}"

def find_mutations(dna_sequence):
    dna_sequence = dna_sequence.strip()
    results_text = ""

    try:
        print(f"Submitting sequence to BLAST: {dna_sequence[:50]}...")  # Show only first 50 chars for brevity
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, word_size=7)
        blast_result = result_handle.read()
        result_handle.close()

        result_file = io.StringIO(blast_result)
        blast_record = NCBIXML.read(result_file)
    except Exception as e:
        return f"An error occurred during BLAST search: {e}\n{traceback.format_exc()}"

    if not blast_record.alignments:
        results_text = "No alignments found. The sequence may not match any known reference."
    else:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 1.0:
                    results_text += "\n**** Alignment Found ****\n"
                    results_text += f"Sequence: {alignment.title}\n"
                    results_text += f"Length: {alignment.length}\n"
                    results_text += f"e-value: {hsp.expect}\n"
                    results_text += f"Identities: {hsp.identities}/{hsp.align_length}\n"
                    results_text += f"Query: {hsp.query}\n"
                    results_text += f"Match: {hsp.match}\n"
                    results_text += f"Subject: {hsp.sbjct}\n"

                    mismatches = [i for i in range(len(hsp.match)) if hsp.match[i] != "|"]
                    if not mismatches:
                        results_text += "No mismatches (no mutations found).\n"
                    else:
                        for idx in mismatches:
                            results_text += f"Position {idx}: Query = {hsp.query[idx]} | Subject = {hsp.sbjct[idx]}\n"

                        for idx in mismatches:
                            query_base = hsp.query[idx]
                            subject_base = hsp.sbjct[idx]
                            position = idx + 1
                            mutation = f"{query_base}{position}{subject_base}"
                            results_text += f"Mutation: {query_base} -> {subject_base} at position {position}\n"

                            disorder_info = search_disorder(mutation)
                            results_text += f"Disorder Info:\n{disorder_info}\n"

    return results_text


class MutationFinderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Mutation Finder")

        self.create_widgets()

    def create_widgets(self):
        self.label = tk.Label(self.root, text="Enter DNA Sequence:")
        self.label.pack(pady=5)

        self.text_area = scrolledtext.ScrolledText(self.root, wrap=tk.WORD, height=10, width=50)
        self.text_area.pack(pady=5)

        self.search_button = tk.Button(self.root, text="Find Mutations", command=self.find_mutations)
        self.search_button.pack(pady=5)

        self.result_area = scrolledtext.ScrolledText(self.root, wrap=tk.WORD, height=15, width=80)
        self.result_area.pack(pady=5)

    def find_mutations(self):
        dna_sequence = self.text_area.get("1.0", tk.END).strip()
        if not dna_sequence:
            messagebox.showwarning("Input Error", "Please enter a DNA sequence.")
            return

        print("Processing sequence...")
        result_text = find_mutations(dna_sequence)
        self.result_area.delete("1.0", tk.END)
        self.result_area.insert(tk.END, result_text)
        print("Finished processing.")

if __name__ == "__main__":
    root = tk.Tk()
    app = MutationFinderApp(root)
    root.mainloop()
