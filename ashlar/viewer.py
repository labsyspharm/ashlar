import napari
import numpy as np
import networkx as nx


_colormaps = ['red', 'green', 'blue']


def view_edges(
        aligner,
        viewer=None,
        tiles=None,
        vmin=None,
        vmax=None,
        gamma=1/2.2,
):
    """View EdgeAligner results in a Napari viewer using a checkerboard style.

    This "checkerboard" visualization colors tiles in pure red, green, and blue
    and applies additive blending to make good or bad alignment easy to spot at
    a per-pixel level. Pure yellow, cyan, and magenta objects indicate good
    alignment, whereas red, green, and blue fringes indicate misalignment.

    All tiles' contrast limits and gamma are linked so you can adjust those
    settings for all tiles by adjusting any one.

    Parameters
    ----------
    aligner: EdgeAligner
        The EdgeAligner to visualize.
    viewer: optional napari.Viewer
        If present, add layers to this viewer instead of creating a new one.
    tiles: optional List[int]
        If present, visualize only the given tiles instead of every one.
    vmin, vmax: optional numeric
        Initial contrast limits. Detected from data if not provided.
    gamma: optional numeric
        Initial gamma correction. Default is appropriate for linear data.

    Returns
    -------
    napari.viewer.Viewer

    """

    if tiles is None:
        tiles = range(aligner.metadata.num_images)
    else:
        if not set(tiles) <= set(range(aligner.metadata.num_images)):
            raise ValueError("tiles contains invalid tile index values")

    node_colors = nx.greedy_color(aligner.neighbors_graph)
    num_colors = max(node_colors.values()) + 1
    if num_colors > 3:
        raise ValueError(f"neighbors_graph requires more than 3 colors")

    dtype = aligner.metadata.pixel_dtype
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
    elif np.issubdtype(dtype, np.floating):
        info = np.iinfo(dtype)
    else:
        raise ValueError(f"Can't display {dtype} image data")
    compute_vmin = vmin is None
    compute_vmax = vmax is None
    if compute_vmin:
        dmin = info.min
        vmin = np.inf
    else:
        dmin = vmin
    if compute_vmax:
        dmax = info.max
        vmax = -np.inf
    else:
        dmax = vmax

    # We wait until now, after we've passed all the argument validations, to
    # create the viewer window. It seems like bad form to throw up an empty
    # window only to then error out due to a bad arg value.
    if viewer is None:
        viewer = napari.Viewer()

    images = []

    # Link contrast sliders for all tiles.
    def contrast_limits_callback(event):
        source_image = event.source
        new_clims = source_image.contrast_limits
        for image in images:
            # this is critical to avoid a RecursionError
            if image._contrast_limits != new_clims:
                image.contrast_limits = new_clims

    # Link gamma sliders for all tiles.
    def gamma_callback(event):
        source_image = event.source
        new_gamma = source_image.gamma
        for image in images:
            # this is critical to avoid a RecursionError
            if image._gamma != new_gamma:
                image.gamma = new_gamma

    for i in tiles:
        data = aligner.reader.read(i, aligner.channel)
        image = viewer.add_image(
            data,
            name=str(i),
            translate=aligner.positions[i],
            colormap=_colormaps[node_colors[i]],
            contrast_limits=(dmin, dmax),
            gamma=gamma,
            blending='additive',
            interpolation='bilinear',
        )
        image.events.contrast_limits.connect(contrast_limits_callback)
        image.events.gamma.connect(gamma_callback)
        images.append(image)
        if compute_vmin:
            vmin = min(vmin, np.min(data))
        if compute_vmax:
            vmax = max(vmax, np.max(data))
    if compute_vmin or compute_vmax:
        for image in images:
            image.contrast_limits = (vmin, vmax)

    bounds_min = aligner.positions[tiles].min(axis=0)
    bounds_max = aligner.positions[tiles].max(axis=0) + aligner.metadata.size
    rect = tuple(bounds_min[::-1]) + tuple(bounds_max[::-1])
    viewer.events.reset_view(rect=rect)

    return viewer
