variable "EXCEPTION_HANDLING" {
  default = "-fwasm-exceptions"
}

variable "BOOST_MAJOR_VERSION" {
  default = "1"
}

variable "BOOST_MINOR_VERSION" {
  default = "87"
}

variable "BOOST_PATCH_VERSION" {
  default = "0"
}

variable "FREETYPE_VERSION" {
  default = "2.13.3"
}

variable "ZLIB_VERSION" {
  default = "1.3.2"
}

variable "EMSDK_VERSION" {
  default = "latest"
}

target "deps" {
  dockerfile = "Dockerfile_1_deps"
  context    = "Code/MinimalLib/docker"
  args = {
    EXCEPTION_HANDLING  = EXCEPTION_HANDLING
    BOOST_MAJOR_VERSION = BOOST_MAJOR_VERSION
    BOOST_MINOR_VERSION = BOOST_MINOR_VERSION
    BOOST_PATCH_VERSION = BOOST_PATCH_VERSION
    FREETYPE_VERSION    = FREETYPE_VERSION
    ZLIB_VERSION        = ZLIB_VERSION
    EMSDK_VERSION       = EMSDK_VERSION
  }
  tags = ["rdkit-minimallib-deps:latest"]
}

target "src" {
  dockerfile = "Code/MinimalLib/docker/Dockerfile_2_rdkit_copy_from_local"
  context    = "."
  contexts = {
    rdkit-minimallib-deps = "target:deps"
  }
  tags = ["rdkit-minimallib-rdkit-src:latest"]
}

target "build" {
  dockerfile = "Code/MinimalLib/docker/Dockerfile_3_rdkit_build"
  context    = "Code/MinimalLib/docker"
  contexts = {
    rdkit-minimallib-rdkit-src = "target:src"
  }
  args = {
    EXCEPTION_HANDLING  = EXCEPTION_HANDLING
    BOOST_MAJOR_VERSION = BOOST_MAJOR_VERSION
    BOOST_MINOR_VERSION = BOOST_MINOR_VERSION
    BOOST_PATCH_VERSION = BOOST_PATCH_VERSION
  }
  tags = ["rdkit-minimallib:latest"]
}

target "export" {
  dockerfile = "Code/MinimalLib/docker/Dockerfile_4_rdkit_export"
  context    = "Code/MinimalLib/docker"
  contexts = {
    rdkit-minimallib = "target:build"
  }
  output = ["type=local,dest=dist"]
}
