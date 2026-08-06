from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field, field_validator, model_validator
from importlib.resources import files

DEFAULTS = {

    "project": "Project",

    "threads": 16,

    "primer":
        "CAAGCGTTGGCTTCTCGCATCT",

    "linker1":
        "GTGGCCGATGTTTCGCATCGGCGTACGACT",

    "linker2":
        "ATCCACGTGCTTGAGAGGCCAGAGCATTCG",
    
    "barcode_file": str(files("most.barcode") / "20240614_2500barcode_AB_update.txt"),

    "rna_barcode_file":
        str(files("most.barcode") / "20240614_2500barcode_AB_update_RNA.txt"),

    "hdist": 3,

    "rna_lib": "illumina",

    "adapter": "AAGCAGTGGTATCAACGCAGAGTGAATGGG",

    "adapter_mismatch": 1
}

def empty_to_default(v, key):

    """
    Replace None / empty string with default value
    """

    if key not in DEFAULTS:
        return v

    if v is None:
        return DEFAULTS[key]

    if isinstance(v, str):

        if v.strip() == "":
            return DEFAULTS[key]

    return v

class SequenceFile(BaseModel):

    model_config = {
        "extra": "forbid"
    }

    file1: str

    file2: str


class Reference(BaseModel):

    model_config = {
        "extra": "forbid"
    }

    genome: str

    chromap_index: Optional[str] = None

    fa_file: str

    gtf_file: str

    star_index: str

    barcode_file: Optional[str] = None

    @field_validator("barcode_file", mode="before")
    @classmethod
    def replace_defaults(cls, v):

        return v

    @field_validator("genome", mode="before")
    @classmethod
    def check_genome(cls, v):
        if v is None:
            raise ValueError("reference.genome is required")

        if isinstance(v, str):
            v = v.strip()
            if not v:
                raise ValueError("reference.genome cannot be empty")
        return v
    
class Advanced(BaseModel):

    model_config = {
        "extra": "forbid"
    }

    primer: str = \
        DEFAULTS["primer"]

    linker1: str = \
        DEFAULTS["linker1"]

    linker2: str = \
        DEFAULTS["linker2"]

    UMI: Optional[str] = None

    BC: Optional[str] = None

    hdist: int = \
        DEFAULTS["hdist"]

    rna_lib: str = \
        DEFAULTS["rna_lib"]

    adapter: str = \
        DEFAULTS["adapter"]

    adapter_mismatch: int = \
        DEFAULTS["adapter_mismatch"]

    @field_validator("*", mode="before")
    @classmethod
    def replace_defaults(cls, v, info):

        return empty_to_default(
            v,
            info.field_name
        )


class Runtime(BaseModel):

    k1: Optional[int] = None

    k2: Optional[int] = None

    bc2_start: Optional[int] = None

    bc2_end: Optional[int] = None

    bc1_start: Optional[int] = None

    bc1_end: Optional[int] = None

    restrictleft1: Optional[int] = None

    restrictleft2: Optional[int] = None

    seq_start: Optional[int] = None

    umi_start: Optional[int] = None

    umi_len: Optional[int] = None


class Config(BaseModel):

    model_config = {
        "extra": "forbid"
    }

    project: str = DEFAULTS["project"]

    method: str

    sequence_file: SequenceFile

    reference: Reference

    out_dir: Optional[str] = None

    threads: int = DEFAULTS["threads"]

    advanced: Advanced = Field(
        default_factory=Advanced
    )

    runtime: Runtime = Field(
        default_factory=Runtime
    )

    tools: dict = Field(
        default_factory=dict
    )

    @field_validator(
        "project",
        "threads",
        mode="before"
    )
    @classmethod
    def replace_defaults(cls, v, info):

        return empty_to_default(
            v,
            info.field_name
        )
    
    @field_validator("out_dir", mode="before")
    @classmethod
    def default_outdir(cls, v, info):

        if v is not None:

            if isinstance(v, str):

                if v.strip() != "":
                    return v

        data = info.data

        if "sequence_file" not in data:
            return "output"

        file1 = data["sequence_file"].file1

        parent = Path(file1).parent
        
        return str(parent / "output")
    
    @field_validator("method", mode="before")
    @classmethod
    def normalize_method(cls, v):

        if v is None:

            raise ValueError(
                "method is required"
            )

        v = v.upper()

        if v not in {
            "RNA",
            "ATAC",
            "ZUMIS",
            "ASTRO",
            "DMT"
        }:
            raise ValueError(
                f"Unknown method: {v}"
            )

        return v
    
    @field_validator("threads")
    @classmethod
    def check_threads(cls, v):

        if v <= 0:
            return DEFAULTS["threads"]

        return v

    @model_validator(mode="after")
    def check_tools(self):
        if not self.reference.barcode_file:
            if self.method == "RNA":
                self.reference.barcode_file = DEFAULTS["rna_barcode_file"]
            else:
                self.reference.barcode_file = DEFAULTS["barcode_file"]

        if self.method == "ZUMIS":
            if not self.tools.get("zumis"):
                raise ValueError(
                    "method=ZUMIS requires tools.zUMIs configuration"
                )

        if self.method == "ASTRO":
            if not self.tools.get("astro"):
                raise ValueError(
                    "method=ASTRO requires tools.ASTRO configuration"
                )

        return self