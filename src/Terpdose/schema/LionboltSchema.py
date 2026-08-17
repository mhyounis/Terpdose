# This file can be edited directly, however, it should always be based upon $LIONBOLT/src/IO/schema.yaml and $NITTANY/src/IO/schema.yaml
# including taking care to keep the schema_version in agreement

SCHEMA_VERSION = '1.0'

ATTR_PROBTYPE = 'problem_type'

GROUP_XSLIBRARY = 'XS_library'

ATTR_OBJTYPE = 'object_type'

GROUP_MESH           = 'mesh'
DATASET_R            = 'global_nodes'
DATASET_CONNECTIVITY = 'connectivity'
DATASET_OFFSET       = 'offset'
DATASET_VOLS         = 'element_volumes'
ATTR_NUMMATS_M       = 'number_of_materials'
ATTR_NUMELS          = 'number_of_elements'
ATTR_NUMSD           = 'number_of_spatial_dofs'
ATTR_NUMKG           = 'number_of_global_nodes'
GROUP_MAT2SD         = 'material_to_spatial_dof'

ATTR_NUMANGLES     = 'number_of_angular_dofs'
DATASET_ANGWEIGHTS = 'angular_weights'
DATASET_ABSCISSAE  = 'abscissae'

ATTR_NUMENERGIES = 'number_of_energy_dofs'
DATASET_EGRID    = 'energy_grid'

IS_PARTICLE       = 'particle'
GROUP_FL          = 'fluence'
GROUP_FL_UNC      = 'uncollided_fluence'
GROUP_ANG_FL      = 'angular_fluence'
GROUP_ANG_FL_UNC  = 'uncollided_angular_fluence'

PATTERN_ENERGY    = 'energy_'
PATTERN_ANGLE     = 'angle_'
PATTERN_MATERIALS = 'material_'

PATTERN_DEP  = '_deposition'
DATASET_EDEP = 'energy'
DATASET_DDEP = 'dose'
DATASET_CDEP = 'charge'
ATTR_EDEP    = 'total_energy'
ATTR_DDEP    = 'total_dose'
ATTR_CDEP    = 'total_charge'