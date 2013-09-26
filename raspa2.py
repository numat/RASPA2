"""
WRAPSA2, the python wrapper for RASPA2.

This is intended 1. to allow users to automate + simplify their workflows and
2. to enable scaling of simulations to millions of structures.
"""
from ctypes import cdll, c_void_p, c_char_p, c_bool, cast
import os
env = os.environ

from RASPA.parser import parse

raspa = cdll.LoadLibrary("/usr/lib/libraspa2.so")
raspa.run.argtypes = (c_char_p, c_char_p, c_char_p, c_bool)
raspa.run.restype = c_void_p


def run(structure, molecule_name, temperature, pressure, helium_void_fraction,
        unit_cells=(1, 1, 1), framework_name="streamed",
        simulation_type="MonteCarlo", cycles=2000, init_cycles=1000,
        forcefield="CrystalGenerator", input_file_type="cif"):
    """Runs a simulation with the specified parameters.

    Args:
        structure: The structure to test for adsorption, as a string of type
            `input_file_type` (default is "cif").
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: The temperature of the simulation, in Kelvin.
        pressure: The pressure of the simulation, in Pascals.
        helium_void_fraction: The helium void fraction of the input structure.
            Can be estimated by running RASPA with helium.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        framework_name: (Optional) If not streaming, this will load the
            structure at `$RASPA_DIR/share/raspa/structures`. Ignored if
            streaming.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """
    return parse(run_script(create_script(**locals()), structure))


def run_script(input_script, structure=None, raspa_dir="auto"):
    """Runs RASPA on the inputted structure, returning simulation results.

    Args:
        input_script: A RASPA simulation script, as a loaded string.
        structure: (Optional) Data encoding a CIF, MOL, or CSSR file. If not
            specified, the program will look for the file specified by
            "FrameworkName" in `$RASPA_DIR/share/raspa/structures`.
        raspa_dir: (Optional) The root directory for RASPA. If unspecified,
            uses the RASPA_DIR environment variable if exists, HOME if not.
    Returns:
        A string representing the output data. (TODO This should be parsed
        and returned as a mutable data structure)

    If a structure is not specified, then nothing will be streamed and this
    will output files and folders. Otherwise, this function will return the
    output of RASPA, as a string.
    """
    if raspa_dir == "auto":
        raspa_dir = env["RASPA2_DIR" if "RASPA2_DIR" in env else "HOME"]
    ptr = raspa.run(input_script, structure or "", raspa_dir, True)
    return cast(ptr, c_char_p).value[:]


def create_script(molecule_name, temperature, pressure, helium_void_fraction,
                  unit_cells=(1, 1, 1), framework_name="streamed",
                  simulation_type="MonteCarlo", cycles=2000,
                  init_cycles=1000, forcefield="CrystalGenerator",
                  input_file_type="cif", **kwargs):
    """Creates a RASPA simulation input file from parameters.

    Args:
        original_cif_file: cif file is read to determine correct unit cell size
            for the simulation
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: The temperature of the simulation, in Kelvin.
        pressure: The pressure of the simulation, in Pascals.
        helium_void_fraction: The helium void fraction of the input structure.
            Can be estimated by running RASPA with helium.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        framework_name: (Optional) If not streaming, this will load the
            structure at `$RASPA_DIR/share/raspa/structures`. Ignored if
            streaming.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
        charged: (Optional) A boolean indicating whether or not to use Ewald
            charge parameters.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """

    he, mol, ift = helium_void_fraction, molecule_name, input_file_type
    return "\n".join(["SimulationType                %s" % simulation_type,
                      "NumberOfCycles                %d" % cycles,
                      "NumberOfInitializationCycles  %d" % init_cycles,
                      "PrintEvery                    %i" % (cycles / 10),
                      "RestartFile                   no",
                      "",
                      "Forcefield                    %s" % forcefield,
                      "CutOff                      12.8",
                      "ChargeMethod                  Ewald",
                      "EwaldPrecision                1e-6",
                      "UseChargesFrom%sFile          yes" % ift.upper(),
                      "",
                      "Framework                     0",
                      "FrameworkName                 %s" % framework_name,
                      "InputFileType                 %s" % ift,
                      "UnitCells                     %d %d %d" % unit_cells,
                      "HeliumVoidFraction            %.2f" % he,
                      "ExternalTemperature           %.1f" % temperature,
                      "ExternalPressure              %d" % pressure,
                      "Movies                        no",
                      "WriteMoviesEvery              100",
                      "",
                      "Component 0 MoleculeName             %s" % mol,
                      "            StartingBead             0",
                      "            MoleculeDefinition       TraPPE",
                      "            IdealGasRosenbluthWeight 1.0",
                      "            TranslationProbability   1.0",
                      "            RotationProbability      1.0",
                      "            ReinsertionProbability   1.0",
                      "            SwapProbability          1.0",
                      "            CreateNumberOfMolecules  0"])


def run_mixture(structure, molecules, mol_fractions, temperature,
                pressure, helium_void_fraction,
                unit_cells=(1, 1, 1), framework_name="streamed",
                simulation_type="MonteCarlo", cycles=2000,
                init_cycles=1000, forcefield="CrystalGenerator",
                input_file_type="cif"):
    """Runs a simulation with mixture of gases.

    Args:
        structure: The structure to test for adsorption, as a string of type
            `input_file_type` (default is "cif").
        molecules: The molecules to test for adsorption. Files of the same
            names must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        mol_fractions: The mol fractions of each gas that you want to separate.
            Corresponds to the `molecules` list.
        temperature: The temperature of the simulation, in Kelvin.
        pressure: The pressure of the simulation, in Pascals.
        helium_void_fraction: The helium void fraction of the input structure.
            Can be estimated by running RASPA with helium.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        framework_name: (Optional) If not streaming, this will load the
            structure at `$RASPA_DIR/share/raspa/structures`. Ignored if
            streaming.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """
    he, ift = helium_void_fraction, input_file_type
    script = "\n".join(["SimulationType                %s" % simulation_type,
                        "NumberOfCycles                %d" % cycles,
                        "NumberOfInitializationCycles  %d" % init_cycles,
                        "PrintEvery                    %i" % (cycles / 10),
                        "RestartFile                   no",
                        "",
                        "Forcefield                    %s" % forcefield,
                        "CutOff                        12.8",
                        "ChargeMethod                  Ewald",
                        "EwaldPrecision                1e-6",
                        "UseChargesFrom%sFile          yes" % ift.upper(),
                        "",
                        "Framework                     0",
                        "FrameworkName                 %s" % framework_name,
                        "InputFileType                 %s" % input_file_type,
                        "UnitCells                     %d %d %d" % unit_cells,
                        "HeliumVoidFraction            %.2f" % he,
                        "ExternalTemperature           %.1f" % temperature,
                        "ExternalPressure              %d" % pressure,
                        "Movies                        no",
                        "WriteMoviesEvery              100",
                        ""])

    changes_list = " ".join(str(n) for n in range(len(molecules)))
    for i, item in enumerate(molecules):
        script += "\n".join(["",
                            "Component %i MoleculeName                 %s"
                            % (i, item),
                            "             StartingBead                 0",
                            "             MoleculeDefinition           TraPPE",
                            "             MolFraction                  %.2f"
                            % mol_fractions[i],
                            "             IdentityChangeProbability    1.0",
                            "                 NumberOfIdentityChanges       %i"
                            % len(molecules),
                            "                 IdentityChangesList           %s"
                            % changes_list,
                            "             IdealGasRosenbluthWeight     1.0",
                            "             TranslationProbability       1.0",
                            "             RotationProbability          1.0",
                            "             ReinsertionProbability       1.0",
                            "             SwapProbability              1.0",
                            "             CreateNumberOfMolecules      0"])
    return parse(run_script(script, structure))


def get_geometric_surface_area(structure, unit_cells=(1, 1, 1), cycles=500,
                               input_file_type="cif", units="m^2/g",
                               forcefield="CrystalGenerator"):
    """Calculates the geometric surface area of an inputted structure.

    Args:
        structure: The structure to use, as a string of type
            `input_file_type` (default is "cif").
        input_file_type: (Optional) The type of input structure. Assumes cif.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        cycles: (Optional) The number of simulation cycles to run.
        units: (Optional) The units in which to return the surface area. Can be
            "m^2/g", "A^2", or "m^2/cm^3".
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The geometric surface area, as a float.
    """
    script = "\n".join(["SimulationType              MonteCarlo",
                        "NumberOfCycles              %d" % cycles,
                        "PrintEvery                  %d" % (cycles / 10),
                        "PrintPropertiesEvery        %d" % (cycles / 10),
                        "",
                        "Forcefield                  %s" % forcefield,
                        "CutOff                      12.8",
                        "",
                        "Framework                   0",
                        "FrameworkName               streamed",
                        "InputFileType               %s" % input_file_type,
                        "UnitCells                   %d %d %d" % unit_cells,
                        "SurfaceAreaProbeDistance    Sigma",
                        "",
                        "Component 0 MoleculeName             N2",
                        "            StartingBead             0",
                        "            MoleculeDefinition       TraPPE",
                        "            SurfaceAreaProbability   1.0",
                        "            CreateNumberOfMolecules  0"])
    output = run_script(script, structure)
    info = parse(output)
    return info["Average Surface Area"]["[%s]" % units][0]


def get_helium_void_fraction(structure, unit_cells=(1, 1, 1), cycles=2000,
                             input_file_type="cif",
                             forcefield="CrystalGenerator"):
    """Calculates the helium void fraction of the inputted structure.

    Args:
        structure: The structure to test for helium void fraction,
            as a string of type 'input_file_type` (default is "cif").
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        cycles: (Optional) The number of simulation cycles to run.
        input_file_type: (Optional) The type of input structure. Assumes cif.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The helium void fraction of the structure, as a float.
    """
    script = "\n".join(["SimulationType              MonteCarlo",
                        "NumberOfCycles              %d" % cycles,
                        "PrintEvery                  %i" % (cycles / 10),
                        "PrintPropertiesEvery        %i" % (cycles / 10),
                        "",
                        "CutOff                      12.8",
                        "Forcefield                  %s" % forcefield,
                        "",
                        "Framework                   0",
                        "FrameworkName               streamed",
                        "InputFileType               %s" % input_file_type,
                        "UnitCells                   %d %d %d" % unit_cells,
                        "ExternalTemperature         298.0",
                        "",
                        "Component 0 MoleculeName             helium",
                        "            MoleculeDefinition       TraPPE",
                        "            WidomProbability         1.0",
                        "            CreateNumberOfMolecules  0"])
    output = run_script(script, structure)
    info = parse(output)
    return info["Average Widom Rosenbluth factor"]["Widom"][0]


def get_density(molecule_name, temperature, pressure, cycles=5000,
                init_cycles=2500, forcefield="CrystalGenerator"):
    """Calculates the density through an NPT ensemble.

    Args:
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: The temperature of the simulation, in Kelvin.
        pressure: The pressure of the simulation, in Pascals.
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The density, as a float, in kg/m^3.
    """
    mol = molecule_name
    script = "\n".join(["SimulationType                   MonteCarlo",
                        "NumberOfCycles                   %d" % cycles,
                        "NumberOfInitializationCycles     %d" % init_cycles,
                        "PrintEvery                       %d" % (cycles / 10),
                        "",
                        "Forcefield                       %s" % forcefield,
                        "",
                        "Box                              0",
                        "BoxLengths                       30 30 30",
                        "ExternalTemperature              %.1f" % temperature,
                        "ExternalPressure                 %d" % pressure,
                        "",
                        "VolumeChangeProbability          0.25",
                        "",
                        "Component 0 MoleculeName               %s" % mol,
                        "            MoleculeDefinition         TraPPE",
                        "            TranslationProbability     0.5",
                        "            ReinsertionProbability     0.5",
                        "            CreateNumberOfMolecules    256"])
    info = parse(run_script(script))
    return info["Average Density"]["[kg/m^3]"][0]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Molecular-dynamic simulator "
                                     "for crystal structures.")
    parser.add_argument("input", type=str, help="A simulation script. Can "
                        "be either a filepath or the input data itself")
    parser.add_argument("--structure", type=str, default="", help="An input "
                        "structure. Can be either a filepath or the input "
                        "data itself")
    args = parser.parse_args()

    # Support both filepaths and data-streamed strings
    try:
        with open(args.input) as in_file:
            args.input = in_file.read()
    except IOError:
        pass
    try:
        with open(args.structure) as in_file:
            args.structure = in_file.read()
    except IOError:
        pass

    print(run_script(args.input, args.structure))
