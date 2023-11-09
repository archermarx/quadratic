project(
    'quadratic', 
    'cpp',
    version: 'v0.1.0',
    default_options: ['cpp_std=c++17']
)

inc = include_directories('include')

subdir('test')

shared_library(
    meson.project_name(),
    [],
    install: true,
    gnu_symbol_visibility: 'hidden',
    include_directories: inc
)

test_exe = executable(
    'runtests', test_sources,
    include_directories: inc,
    install: false
)

test('tests', test_exe)