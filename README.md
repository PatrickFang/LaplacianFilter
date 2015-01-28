- How to compile the code: simply use make command under the program directory.

- How to run the code:       the program can take up to 4 arguments in the following order
                                     -- image path, i.e: while under the program directory, "input_png/flower.png"
                                     -- sigma_r (double value)
                                     -- alpha     (double value)
                                     -- beta       (double value)
                                    However, for simplicity in testing, I have coded the program to give default values to the last three double values while they are not input.


- DISPLAY_WINDOW preprocessor option, locating at line 23.

- Can use the Preprocessor option of resizing the input image to its 1/4 size, which locates at line 24 of main. to reduce run time

- The program only takes coloured image, so the ones with RBG, 3 channels
