from matplotlib import pyplot as plt
import subprocess, os

def plot_as_emf(figure, **kwargs):
    inkscape_path = kwargs.get('inkscape', "C://Program Files//Inkscape//inkscape.exe")
    filepath = kwargs.get('filename', None)

    if filepath is not None:
        path, filename = os.path.split(filepath)
        filename, extension = os.path.splitext(filename)

        svg_filepath = os.path.join(path, filename+'.svg')
        emf_filepath = os.path.join(path, filename+'.emf')

        figure.savefig(svg_filepath, format='svg')

        subprocess.call([inkscape_path, svg_filepath, '--export-emf', emf_filepath])
        os.remove(svg_filepath)



