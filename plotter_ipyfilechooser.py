from ipyfilechooser import FileChooser
import os

# Create and display a FileChooser widget
fc = FileChooser()
fc.default_path = os.path.expanduser('~')
fc.use_dir_icons = True
fc.title = '<b>Select File</b>'
display(fc)

# Print the selected path, filename, or both
print(fc.selected_path)
print(fc.selected_filename)
print(fc.selected)

# Change defaults and reset the dialog
fc.default_path = '/Users/crahan/'
fc.default_filename = 'output.txt'
fc.reset()

# Shorthand reset
fc.reset(path='/Users/crahan/', filename='output.txt')

# Change hidden files
fc.show_hidden = True

# Show or hide folder icons
fc.use_dir_icons = True

# Change the title (use '' to hide)
fc.title = '<b>FileChooser title</b>'

# Sample callback function
def change_title(chooser):
    chooser.title = '<b>Callback function executed</b>'

# Register callback function
fc.register_callback(change_title)