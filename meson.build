project('c++ tutorial', 'cpp')

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