import os
import PIL

from distutils.spawn import find_executable

meme_path = "/%s" % find_executable('meme').strip("/meme")

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def mkdir_p(path):
    import os
    import errno

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tab2fasta(posneg, fasta, background):
    i = 0
    # print fasta
    with open(posneg) as tab:
        with open(fasta, 'w') as fa:
            with open(background, 'w') as bg:
                for line in tab:
                    details = line.split()
                    if len(details) == 2:
                        pos = 1
                    else:
                        pos = 2
                    if i < 500:
                        fa.write(">" + line.split()[0] + '\n' + line.split()[pos] + "\n")
                    else:
                        bg.write(">" + line.split()[0] + '\n' + line.split()[pos] + "\n")
                    i += 1


def rotate_image(figure_in, figure_out):
    src_im = PIL.Image.open(figure_in)
    im = src_im.rotate(270, expand=True)
    im.save(figure_out)


revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in x][::-1])