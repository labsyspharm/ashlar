import collections
import datetime
import json
import pathlib
import re
from dataclasses import dataclass

import lxml.etree
import ome_types


def read_metadata_file(metadata_path=None, rcjob_path=None):
    metadata = {}
    rcjob = {}

    if metadata_path is not None:
        with open(metadata_path) as f:
            metadata = f.read().strip().split("\n")
        metadata = {ii.split("=")[0]: ii.split("=")[1] for ii in metadata if "=" in ii}

    if rcjob_path is not None:
        with open(rcjob_path) as f:
            rcjob = json.load(f)

    metadata = collections.defaultdict(lambda: None, metadata)
    rcjob = collections.defaultdict(lambda: None, rcjob)

    mag = metadata["BhInstrLensName"]
    if mag is not None:
        mag = float(mag.lower().replace("x", ""))

    objective = ome_types.model.Objective(
        serial_number=metadata["BhInstrLensId"],
        lens_na=metadata["BhInstrLensNA"],
        nominal_magnification=mag,
    )
    microscope = ome_types.model.Microscope(
        manufacturer="RareCyte",
        model="CyteFinder II HT",
        serial_number=rcjob["jobHostname"],
    )

    instrument = ome_types.model.Instrument(
        microscope=microscope, objectives=[objective]
    )

    job_datetime = rcjob["jobSubmissionTime"]
    if job_datetime is not None:
        job_datetime = datetime.datetime.strptime(job_datetime, "%Y%m%d_%H%M%S_%f")

    n_channels = len([kk for kk in metadata.keys() if re.match("Biomarker\d+Em", kk)])
    channels = []
    planes = []
    for i in range(n_channels):
        i += 1
        channels.append(
            ome_types.model.Channel(
                id=f"Channel:{i-1}",
                emission_wavelength=metadata[f"Biomarker{i}Em"],
                excitation_wavelength=metadata[f"Biomarker{i}Ex"],
                name=metadata[f"Biomarker{i}Name"],
            )
        )
        planes.append(
            ome_types.model.Plane(
                the_c=i - 1,
                the_t=0,
                the_z=0,
                exposure_time=metadata[f"Biomarker{i}Exp"],
                exposure_time_unit="s",
            )
        )
    pixel = ome_types.model.Pixels(
        dimension_order="XYCZT",
        size_c=n_channels,
        size_t=1,
        size_z=1,
        size_x=1,
        size_y=1,
        type="uint16",
        metadata_only=True,
        channels=channels,
        planes=planes,
    )
    image = ome_types.model.Image(pixels=pixel, acquisition_date=job_datetime)

    return ome_types.model.OME(instruments=[instrument], images=[image])


@dataclass
class ChannelMetadata:
    """Class for keeping track of channel metadata"""

    image_filepath: str | pathlib.Path = None
    ashlar_version: str = None
    alignment_filepath: str | pathlib.Path = None
    acquisition_date: datetime.datetime = None
    is_subtracted: bool = None
    subtraction_filepath: str | pathlib.Path = None
    subtraction_alignment_filepath: str | pathlib.Path = None
    subtraction_instrument: ome_types.model.Instrument = None
    subtraction_channel: ome_types.model.Channel = None
    subtraction_exposure_time: float = None

    def to_xml(self):
        root = lxml.etree.Element("ChannelMetadata")
        metadata = self.__dict__
        for kk, vv in metadata.items():
            if vv is None:
                continue
            kk = "".join(tt.capitalize() for tt in kk.split("_"))
            ke = lxml.etree.Element(kk)
            if hasattr(vv, "to_xml"):
                ke_child = vv.to_xml(
                    include_namespace=False, include_schema_location=False
                )
                ke_child = lxml.etree.fromstring(ke_child)
                ke.append(ke_child)
            elif hasattr(vv, "isoformat"):
                ke.text = vv.isoformat()
            else:
                ke.text = f"{vv}"
            root.append(ke)
        for el in root.findall(".//*[@ID]"):
            # el.attrib["ID"] = ""
            el
        return lxml.etree.tostring(root, pretty_print=True, encoding="unicode")


def _assemble_metadata(ref_path, path, output_path, channels, _is_subtracted=False):
    path = pathlib.Path(path).absolute()
    metadata_path, rcjob_path = _load_metadata_paths(path.parent)
    metadata_ome = read_metadata_file(
        metadata_path=metadata_path, rcjob_path=rcjob_path
    )
    img_ome = ome_types.from_tiff(output_path)

    img_ome.instruments = metadata_ome.instruments

    img_ome.images[0].acquisition_date = metadata_ome.images[0].acquisition_date

    ome_channels = metadata_ome.images[0].pixels.channels
    ome_planes = metadata_ome.images[0].pixels.planes

    if channels is not None:
        ome_channels = [metadata_ome.images[0].pixels.channels[ii] for ii in channels]
        ome_planes = [metadata_ome.images[0].pixels.planes[ii] for ii in channels]

    img_ome.images[0].pixels.channels = ome_channels
    img_ome.images[0].pixels.planes = ome_planes

    for idx, pp in enumerate(img_ome.images[0].pixels.planes):
        pp.the_c = idx

    if ref_path is not None:
        ref_path = pathlib.Path(ref_path).absolute()

    channel_metadata = ChannelMetadata(
        image_filepath=path,
        alignment_filepath=ref_path,
        acquisition_date=img_ome.images[0].acquisition_date,
        is_subtracted=_is_subtracted,
    )
    annotation = ome_types.model.XMLAnnotation(value=channel_metadata.to_xml())
    img_ome.structured_annotations.append(annotation)
    for channel in img_ome.images[0].pixels.channels:
        channel.annotation_refs.append(ome_types.model.AnnotationRef(id=annotation.id))
    return img_ome


def _subtract_metadata(
    bg_path,
    bg_ref_path,
    ab_path,
    ab_ref_path,
    output_path,
    fiducial_channel,
    subtraction_config,
):
    ab_path = pathlib.Path(ab_path).absolute()
    bg_path = pathlib.Path(bg_path).absolute()

    if ab_ref_path is not None:
        ab_ref_path = pathlib.Path(ab_ref_path).absolute()
    if bg_ref_path is not None:
        bg_ref_path = pathlib.Path(bg_ref_path).absolute()

    ome = _assemble_metadata(
        ab_ref_path, ab_path, output_path, channels=None, _is_subtracted=None
    )

    metadata_path, rcjob_path = _load_metadata_paths(bg_path.parent)
    metadata_ome = read_metadata_file(
        metadata_path=metadata_path, rcjob_path=rcjob_path
    )
    metadata_pixels = metadata_ome.images[0].pixels
    for channel in metadata_pixels.channels:
        parts = channel.id.split(":")
        parts.insert(1, "background")
        channel.id = ":".join(parts)

    metadata = {
        "fiducial": dict(is_subtracted=False),
        "subtracted": dict(
            is_subtracted=True,
            subtraction_filepath=bg_path,
            subtraction_alignment_filepath=bg_ref_path,
            subtraction_instrument=metadata_ome.instruments[0],
        ),
    }

    annot_fiducial = ome_types.model.XMLAnnotation(
        value=ChannelMetadata(**metadata["fiducial"]).to_xml()
    )
    annot_subtracted = ome_types.model.XMLAnnotation(
        value=ChannelMetadata(**metadata["subtracted"]).to_xml()
    )
    ome.structured_annotations.extend([annot_fiducial, annot_subtracted])

    ab_channel_ids = [ch_ab["channel_index"] for ch_ab in subtraction_config]
    if fiducial_channel is not None:
        assert fiducial_channel in ab_channel_ids
        subtraction_config[ab_channel_ids.index(fiducial_channel)]["bg_channel"] = {}

    for channel, config in zip(ome.images[0].pixels.channels, subtraction_config):
        bg_channel_id = config["bg_channel"].get("channel_index", None)
        if bg_channel_id is None:
            channel.annotation_refs.append(
                ome_types.model.AnnotationRef(id=annot_fiducial.id)
            )
            continue
        annot_channel = ome_types.model.XMLAnnotation(
            value=ChannelMetadata(
                subtraction_exposure_time=metadata_pixels.planes[
                    bg_channel_id
                ].exposure_time,
                subtraction_channel=metadata_pixels.channels[bg_channel_id],
            ).to_xml()
        )
        ome.structured_annotations.append(annot_channel)
        channel.annotation_refs.extend(
            [
                ome_types.model.AnnotationRef(id=annot_subtracted.id),
                ome_types.model.AnnotationRef(id=annot_channel.id),
            ]
        )
    return ome


def _combine_metadata(input_files, output_path, dna_file_index, dna_channel_number):
    input_omes = [ome_types.from_tiff(pp) for pp in input_files]
    ome = ome_types.from_tiff(output_path)
    ome.instruments = input_omes[0].instruments
    ome.images[0].acquisition_date = input_omes[0].images[0].acquisition_date
    ome.creator = input_omes[0].creator

    ome_channels = []
    ome_planes = []
    ome_structured_annotations = []

    for idx, oo in enumerate(input_omes):
        # add ashlar version annotation
        annot_version = ome_types.model.XMLAnnotation(
            value=ChannelMetadata(ashlar_version=oo.creator).to_xml()
        )
        oo.structured_annotations.append(annot_version)
        for channel in oo.images[0].pixels.channels:
            channel.annotation_refs.append(
                ome_types.model.AnnotationRef(id=annot_version.id)
            )

        # add additional digit to ID numberings
        root = lxml.etree.fromstring(oo.to_xml())
        for el in root.findall(".//*[@ID]"):
            parts = el.attrib["ID"].split(":")
            parts.insert(1, str(idx))
            el.attrib["ID"] = ":".join(parts)
        oo = ome_types.from_xml(lxml.etree.tostring(root))
        _annots = oo.structured_annotations
        channels = oo.images[0].pixels.channels
        planes = oo.images[0].pixels.planes
        if idx != dna_file_index:
            channels.pop(dna_channel_number)
            planes.pop(dna_channel_number)
        ome_channels.extend(channels)
        ome_planes.extend(planes)

        annot_ids = set(
            [
                ann.id
                for channel in oo.images[0].pixels.channels
                for ann in channel.annotation_refs
            ]
        )
        annot = list(filter(lambda x: x.id in annot_ids, _annots))
        ome_structured_annotations.extend(annot)

    for idx, plane in enumerate(ome_planes):
        plane.the_c = idx

    ome.images[0].pixels.channels = ome_channels
    ome.structured_annotations.extend(ome_structured_annotations)
    ome.images[0].pixels.planes = ome_planes
    return ome


def _load_metadata_paths(path):
    rcjob_paths = sorted(path.glob("*.rcjob"))
    if len(rcjob_paths) == 0:
        return None, None
    rcjob_path = rcjob_paths[-1]
    metadata_path = path / rcjob_path.name.replace(".rcjob", ".metadata")
    if not metadata_path.exists():
        metadata_path = None
    return metadata_path, rcjob_path
