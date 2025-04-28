"""
Interactive phylogenetic tree using ete3 and PyQt5
With interactive leaf names: click to print name in console
"""
from PyQt5 import QtCore
from PyQt5.QtWidgets import QGraphicsSimpleTextItem
from ete3 import faces, TreeStyle, NodeStyle, NCBITaxa, TextFace


class TaxoTree():
    """Class that creates an ete3 tree and adds a style to it for later interactive display.

    This class contains the functions necessary to create an ete3 interactive tree and a style
    to it. The tree is created based on the information stored in the controller object and the
    NCBITaxa package wich contains taxonomical information of all species in NCBI database. The
    taxonomical information is then used to build a tree without branch lengths. 
    The style given to the tree is based on the biological name and the percentage of identity.
    Leaves contain an interactive element in form of the identity percentage text that display
    an alignment between all the sequences found for a given organism and the barcode blasted.
    Leaves also have a background colour that depends on the percentage of identity.

    Attributes:
        target (str): name of the organism being studied 
        controller (MainController): controller which will create an instance of the tree
        species_dict (dict): dictionary that contains the name and percentage of identity of each species
        tree (Tree object): tree created by ete3 toolkit
        style (Tree style object): style created by ete3 toolkit
    """
    def __init__(self, species_info, target, controller):
        self.target = target
        self.controller = controller
        self.species_dict = self._get_dict_from_info(species_info)
        self.tree = self._create_tree()
        self.style = self._add_style()

    def _get_dict_from_info(self, species_info):
        """
        Private fucntion to filter the data used to create the tree to only take into account
        species names and percentage of identity.
        """
        species_dict = {}
        for organism in species_info:
            name = organism["Scientific name"]
            id_pct = organism["Identity percentage"]
            species_dict[name] = id_pct
        return species_dict

    def _create_tree(self):
        """
        Private function that creates the ete3 toolkit tree object using the NCBITaxa package to
        obtain the taxonomical information about the species inputted.
        """
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
        # Annotate each node of the tree with the scientific name
        ncbi.annotate_tree(tree, taxid_attr='name')
        # Delete empty nodes without leaves or leaves without identity information
        self._delete_empty_nodes(tree)
        return tree
        
    def _add_style(self):
        """
        Private function that adds a ete3 toolkit style based to the previously generated tree
        """
        # Create style
        ts = TreeStyle()
        # Follow style guidelinges in master_layout()
        ts.layout_fn = self._master_layout
        # Add title with the name of the studied organism
        ts.title.add_face(
            faces.TextFace(self.target, fsize=15),
            column=0
        )
        return ts

    def _master_layout(self, node):
        """
        Private function that defines the style of the tree. 
        If a node is not a leaf it adds a label with the biological name of the rank in blue 
        and italic font to distinguisf from species.
        If a node IS a leaf, it adds the biological name, an interactive element and a bakcground
        colour.
        """
        # Use the scientific name of the node as the general name of the node
        node.name = node.sci_name
        # Check if a node is leaf 
        if not node.is_leaf():
            # Add biological name of the rank in blue and italics
            label = TextFace(node.sci_name, fsize=9, fstyle="italic", fgcolor="blue")
            node.add_face(label, column=0, position="branch-right")
        else:
            # Check if the name exists in the species dictionary to avoid erros with identity percentage
            if node.sci_name in self.species_dict.keys():
                # Get identity percentage from species dictionary
                id_pct = self.species_dict[node.sci_name]
                node.add_feature('id_pct', id_pct)
                # Set the node style to add a background colour depending on the percentage of identity
                node.set_style(NodeStyle())
                node.img_style["bgcolor"] = self._get_colour_by_pct(id_pct, node)
                # Add an interactive text to the leaf so user can click on it 
                face = faces.DynamicItemFace(self._leaf_name_face)
                faces.add_face_to_node(face, node, column=1, position="aligned")         
    
    def _leaf_name_face(self, node):
        """
        Private face generator: returns an InteractiveText item showing the leaf name.
        """
        # Get identity percentage with only two decimal figures
        id_pct = round(self.species_dict[node.name], 2)
        return InteractiveText(f'  {id_pct}', node, self.controller)
    
    def _get_colour_by_pct(self, pct, node):
        """
        Private function that returns a colour in hex format based on the identity percentage. 
        The colour scale goes from red (for 100%) to green (0%) going through yellow. 
        The scale is broken into two sections to present a bigger colour difference between 
        percentages over and under 98%.
        If the leaf contains the name of the target organism, the background colour is blue
        """
        # Check if leaf contains target organism
        if node.name.lower() != self.target.lower():
            blue = 120
            red = 255 if pct >= 98 else int(pct*2.5)
            green = (255-int((pct-98)*125)) if pct >= 98 else 255
        else:
            # Light lue
            red = 100
            green = 215
            blue = 255
        # Convert to hex format and return
        return "#{:02X}{:02X}{:02X}".format(red, green, blue)

    def _delete_empty_nodes(self, tree):
        """
        Deletes leaves whose scientific name does not appear in the species dictionary
        to avoid errors and blank leaves.
        """
        for node in tree.traverse():
            if node.is_leaf() and node.sci_name not in self.species_dict.keys():
                node.delete()
    
    def show_tree(self):
        """
        Creates a resizable window to display the taxonomy tree.
        """
        self.tree.show(tree_style=self.style)


class InteractiveText(QGraphicsSimpleTextItem):
    """Interactive object in text format which is appended to the tree leaves to display more info.
    
    This class inherits from a PyQt5 class to create a custom interactive text that is used to display
    the percentage of identity in a leaf as well as have some properties when hovered and clicked.
    When the text is hovered the text becomes bold to show the user where is the mouse pointing at and
    when it is clicked, the controller is called to create a pop up window with alignments between all
    the sequences under the same species name and the barcode blasted.

    Attributes:
        node (Node object): node where the instance of the class belongs to
        controller (MainController): controller used to call function to create alignment pop up
    """
    def __init__(self, text, node, controller):
        super(InteractiveText, self).__init__(text)
        self.node = node
        self.controller = controller
        # Enable click and hover
        self.setAcceptHoverEvents(True)
        self.setCursor(QtCore.Qt.PointingHandCursor)

    def rect(self):
        """
        Provide rect() for DynamicItemFace compatibility
        """
        return self.boundingRect()

    def mousePressEvent(self, event):
        """
        Event when the text is clicked
        """
        # Instruct controller to show alignment pop up
        self.controller.show_alignment_from_tree(self.node.name)
        super(InteractiveText, self).mousePressEvent(event)

    def hoverEnterEvent(self, event):
        """
        Event when the mouse is hovering the text
        """
        # Make text bold on hover
        font = self.font()
        font.setBold(True)
        self.setFont(font)
        super(InteractiveText, self).hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        """
        Event when the mouse is not hovering the text
        """
        # Revert text style when not hovered
        font = self.font()
        font.setBold(False)
        self.setFont(font)
        super(InteractiveText, self).hoverLeaveEvent(event)

