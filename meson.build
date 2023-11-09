project('quadratic', 'cpp', version: 'v0.1.0', default_options: ['cpp_std=c++20'])

inc = include_directories('include')

subdir('src')
subdir('include')
subdir('test')

executable(
    'main', sources,
    include_directories: inc
)

test_exe = executable(
    'runtests', test_sources,
    include_directories: inc
)

test('tests', test_exe)