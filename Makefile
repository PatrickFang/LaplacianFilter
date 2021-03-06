l: main.o numlevels.o downsample.o child_window.o laplacian_pyramid.o upsample.o reconstruct_laplacian_pyramid.o
	g++ -g main.o numlevels.o downsample.o child_window.o laplacian_pyramid.o upsample.o reconstruct_laplacian_pyramid.o -o l -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
main.o: main.cpp
	g++ -c main.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
numlevels.o: numlevels.cpp numlevels.h 
	g++ -c numlevels.h numlevels.cpp 
downsample.o: downsample.cpp downsample.h 
	g++ -c downsample.cpp downsample.h -O2 -L/usr/X11R6/lib -lm -lpthread -lX11
child_window.o: child_window.cpp child_window.h
	g++ -c child_window.cpp child_window.h
laplacian_pyramid.o: laplacian_pyramid.cpp laplacian_pyramid.h
	g++ -c laplacian_pyramid.cpp laplacian_pyramid.h
upsample.o: upsample.cpp upsample.h
	g++ -c upsample.cpp upsample.h
reconstruct_laplacian_pyramid.o: reconstruct_laplacian_pyramid.cpp reconstruct_laplacian_pyramid.h
	g++ -c reconstruct_laplacian_pyramid.cpp reconstruct_laplacian_pyramid.h
clean:
	$(RM) *.o *~ $(main)
