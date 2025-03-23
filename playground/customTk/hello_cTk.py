import customtkinter as ctk

class BlastButtonFrame(ctk.CTkFrame):
    def __init__(self, parent, blast_params):
        super().__init__(parent)
        self.blast_params = blast_params # store reference to blast checkbox frame

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        # BLAST button
        self.BLAST_button = ctk.CTkButton(self, text='BLAST', command=self.button_callback) # Button
        self.BLAST_button.grid(row=0, column=0, padx=20, pady=5, sticky="nsew") # Grid layout

    def button_callback(self):
        checked_checkboxes = self.blast_params.get_cb()
        print(f'Checked checkboxes: {checked_checkboxes}')



class BlastParametersCheckboxFrame(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent)

        # BLAST parameters checkbox
        self.megablast_checkbox = ctk.CTkCheckBox(self, text='Highly similar sequences (megablast)')
        self.megablast_checkbox.grid(row=0, column=0, padx=(20,5), pady=5, sticky="w")

        self.discontiguous_megablast_checkbox = ctk.CTkCheckBox(self, text='More dissimilar sequences (discontiguous megablast)')
        self.discontiguous_megablast_checkbox.grid(row=1, column=0, padx=(20,5), pady=5, sticky="w")

        self.blastn_checkbox = ctk.CTkCheckBox(self, text='Somewhat similar sequences (blastn)')
        self.blastn_checkbox.grid(row=2, column=0, padx=(20,5), pady=5, sticky="w")

    def get_cb(self):
        checked_checkboxes = []
        if self.megablast_checkbox.get():
            checked_checkboxes.append(self.megablast_checkbox.cget('text'))
        if self.discontiguous_megablast_checkbox.get():
            checked_checkboxes.append(self.discontiguous_megablast_checkbox.cget('text'))
        if self.blastn_checkbox.get():
            checked_checkboxes.append(self.blastn_checkbox.cget('text'))
        return checked_checkboxes
        


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        self.title('BLOOM customTkinter POC')
        self.geometry('400x150')

        # Grid layout row for the whole window
        self.grid_rowconfigure(0, weight=1)

        # Grid layout column for the whole window

        
         # BLAST parameter frame
        self.blast_param_frame = BlastParametersCheckboxFrame(self)
        self.blast_param_frame.grid(row=2, column=0, padx=5, pady=5, sticky="nsw")

        # BLAST button frame
        self.blast_frame = BlastButtonFrame(self, self.blast_param_frame)
        self.blast_frame.grid(row=1, column=0, padx=5, pady=5, sticky="nswe")
        
       

    
    

app = App()
app.mainloop()


