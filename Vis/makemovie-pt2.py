import os

name = 'foo'
framerate = 20 # fps
frames = 150

for frame in range(frames):
    os.system('mv {0}-{1}.tga {0}-{1:04d}.tga'.format(name, frame))

os.system('ffmpeg -r {0} -i {1}-%4d.tga -vb 20M foo2.mpg'.format(framerate,
    name))
