{
  'targets': [
    {
      'target_name': 'chem-trep',
      'sources': [ 'src/napi.chem.cc',
					'src/single_atom.cpp',
					'src/single_bond.cpp',
					'src/molutils.cpp',
					'src/scalar.cpp',
					'src/svg_molecule.cpp',
					'src/mol_struct_common.cpp',
					'src/simple_molecule.cpp',
					'src/edited_molecule.cpp' ],
      'include_dirs': ["<!@(node -p \"require('node-addon-api').include\")"],
      'dependencies': ["<!(node -p \"require('node-addon-api').gyp\")"],
      'cflags!': [ '-fno-exceptions' ],
      'cflags_cc!': [ '-fno-exceptions' ],
      'xcode_settings': {
        'GCC_ENABLE_CPP_EXCEPTIONS': 'YES',
        'CLANG_CXX_LIBRARY': 'libc++',
        'MACOSX_DEPLOYMENT_TARGET': '10.7'
      },
      'msvs_settings': {
        'VCCLCompilerTool': { 'ExceptionHandling': 1 },
      }
    }
  ]
}
