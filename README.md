# Efficient Ceremony

This is a course work for holding efficienct ceremony for Groth16.
Simply do '$cargo run'. It will test the time cost of holding ceremonies of different scales.
The test will last for quite some time since the program will also update the common reference strings in linear time. But the verification increases in logarithmic time.