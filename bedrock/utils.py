import seaborn as sns

def make_palette(categories, palette="colorblind", fallback = "hls"):
    categories = sorted(set(categories))
    n = len(categories)
    
    try:
        colours = sns.color_palette(palette, n)
    except ValueError:
        colours = sns.color_palette(fallback, n)

    return dict(zip(categories, colours))
