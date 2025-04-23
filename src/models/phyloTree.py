"""
Interactive phylogenetic tree using ete3 and PyQt5
With interactive leaf names: click to print name in console
"""
from PyQt5 import QtCore
from PyQt5.QtWidgets import QGraphicsSimpleTextItem
from ete3 import faces, TreeStyle, NodeStyle, NCBITaxa, TextFace


class PhyloTree:
    def __init__(self, species_info, target, controller):
        self.target = target
        self.controller = controller
        self.species_dict = self.get_dict_from_info(species_info)
        self.tree = self.create_tree()
        self.style = self.add_style()

    def get_dict_from_info(self, species_info):
        species_dict = {}
        for organism in species_info:
            name = organism["Scientific name"]
            id_pct = organism["Identity percentage"]
            seq = organism["Subject sequence"]
            species_dict[name] = [id_pct, seq]
        return species_dict

    def create_tree(self):
        # Load phyolgentic info from NCBI
        ncbi = NCBITaxa()
        # Get names of species
        names = list(self.species_dict.keys())
        # Translate names to taxonomic IDs
        name2taxid = ncbi.get_name_translator(names)
        # return a list of lists
        taxid_list = list(name2taxid.values())
        # Flatten the  list of lists into a single list of integers
        taxid_int_list = []
        for element in taxid_list:
            num = int(element[0])
            taxid_int_list.append(num)
        # Get topology returns a tree structure containing the input taxa
        tree = ncbi.get_topology(taxid_int_list)
        ncbi.annotate_tree(tree, taxid_attr='name')
        self.delete_empty_nodes(tree)
        return tree
        
    def add_style(self):
        ts = TreeStyle()
        ts.layout_fn = self.master_layout
        ts.title.add_face(
            faces.TextFace(self.target, fsize=15),
            column=0
        )
        return ts

    def master_layout(self, node):
        """Layout: attach leaf_name_face to each leaf node."""
        node.name = node.sci_name
        # Only add a face label for internal nodes
        if not node.is_leaf():
            label = TextFace(node.sci_name, fsize=9, fstyle="italic", fgcolor="blue")
            node.add_face(label, column=0, position="branch-right")
        else:
            if node.sci_name in self.species_dict.keys():
                id_pct = self.species_dict[node.sci_name][0]
                seq = self.species_dict[node.sci_name][1]
                node.add_feature('id_pct', id_pct)
                node.add_feature('seq', seq) 
                node.set_style(NodeStyle())
                node.img_style["bgcolor"] = self.get_color_by_pct(id_pct, node)
                face = faces.DynamicItemFace(self.leaf_name_face)
                faces.add_face_to_node(face, node, column=1, position="aligned")         
    
    def leaf_name_face(self, node):
        """Face generator: returns an InteractiveText item showing the leaf name."""
        id_pct = round(self.species_dict[node.name][0], 2)
        return InteractiveText(f'  {id_pct}', node, self.controller)

    def show_tree(self):
        self.tree.show(tree_style=self.style)
    
    def get_color_by_pct(self, pct, node):
        if node.name.lower() != self.target.lower():
            blue = 120
            red = 255 if pct >= 98 else int(pct*2.5)
            green = (255-int((pct-98)*125)) if pct >= 98 else 255
        else:
            red = 100
            green = 215
            blue = 255
        return "#{:02X}{:02X}{:02X}".format(red, green, blue)

    def delete_empty_nodes(self, tree):
        for node in tree.traverse():
            if node.is_leaf() and node.sci_name not in self.species_dict.keys():
                node.delete()
                


class InteractiveText(QGraphicsSimpleTextItem):
    def __init__(self, text, node, controller):
        super(InteractiveText, self).__init__(text)
        self.node = node
        self.controller = controller
        # Enable click and hover
        self.setAcceptHoverEvents(True)
        self.setCursor(QtCore.Qt.PointingHandCursor)

    def rect(self):
        # Provide rect() for DynamicItemFace compatibility
        return self.boundingRect()

    def mousePressEvent(self, event):
        # Print leaf name on click
        #self.controller.print_something(self.node.seq)
        self.controller.show_alignment_from_tree(self.node.name)
        super(InteractiveText, self).mousePressEvent(event)

    def hoverEnterEvent(self, event):
        # Optional: make text bold on hover
        font = self.font()
        font.setBold(True)
        self.setFont(font)
        super(InteractiveText, self).hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        # Revert text style when not hovered
        font = self.font()
        font.setBold(False)
        self.setFont(font)
        super(InteractiveText, self).hoverLeaveEvent(event)

