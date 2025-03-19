from ete3 import Tree, TreeStyle, NodeStyle
from Bio import Entrez
import time
import os

# Set your email for NCBI access
Entrez.email = "tiburosodrilos@gmail.com"


def get_taxonomy(species_name):
    try:
        # Search for taxonomy ID
        handle = Entrez.esearch(db="taxonomy", term=species_name)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print(f"Species '{species_name}' not found in NCBI Taxonomy.")
            return None

        tax_id = record['IdList'][0]

        # Get taxonomy lineage
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        lineage = records[0]['Lineage'].split("; ")
        lineage.append(species_name)  # Include species at the end of the lineage

        return lineage
    except Exception as e:
        print(f"Error retrieving taxonomy for {species_name}: {e}")
        return None


def build_tree(species_list):
    tree = Tree()

    for species in species_list:
        lineage = get_taxonomy(species)
        if lineage:
            current_node = tree
            for rank in lineage:
                # Search for existing node, otherwise create a new one
                found = False
                for child in current_node.children:
                    if child.name == rank:
                        current_node = child
                        found = True
                        break
                if not found:
                    new_node = current_node.add_child(name=rank)
                    current_node = new_node
            time.sleep(0.5)  # To avoid NCBI request limits

    return tree


def style_tree(tree):
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "r"
    ts.scale = 0.1

    species_colors = {
        "Homo sapiens": "green",
        "Pan troglodytes": "orange",
        "Mus musculus": "blue",
        "Canis lupus familiaris": "purple",
        "Gallus gallus": "red"
    }

    for node in tree.traverse():
        nstyle = NodeStyle()
        if node.name in species_colors:
            nstyle["fgcolor"] = species_colors[node.name]
        else:
            nstyle["fgcolor"] = "darkblue"

        nstyle["size"] = 1 if node.is_leaf() else 6
        node.set_style(nstyle)

    return ts


def main():
    species_list = [
        "Homo sapiens",
        "Pan troglodytes",
        "Mus musculus",
        "Canis lupus familiaris",
        "Gallus gallus"
    ]

    print("Building phylogenetic tree...")
    tree = build_tree(species_list)

    

    if tree:
        print("\nPhylogenetic Tree:")
        print(tree.get_ascii(show_internal=True))

        ts = style_tree(tree)
        tree.show(tree_style=ts)  # Display the tree graphically
        filename = os.path.expanduser('~/Downloads/tree.png')
        tree.render(filename, tree_style=ts)
        print(f"Tree exported as {filename}")


if __name__ == "__main__":
    main()
