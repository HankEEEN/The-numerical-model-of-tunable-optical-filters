`timescale 1ns/1ns
module my_waveform_gem_tb;


    reg clk;
    reg reset;
    reg signed [8:0] control;
    wire [7:0] my_sine_out;
    wire [7:0] my_cosine_out;
    integer outfile;


    // The instantiation of my waveform generator
    my_waveform_gem uut (
        .clk(clk),
        .reset(reset),
        .control(control),
        .my_sine_out(my_sine_out),
        .my_cosine_out(my_cosine_out)
    );

    

    // Create a file to store the sine and cosine outputs
    initial begin
        outfile = $fopen("output_data.txt", "w");
    end


    // 100 MHz clock - 10 ns clock period (clock cycle)
    initial clk = 0;
    always #5 clk = ~clk;

    initial begin
    // Dump the simulation result to a VCD file
        $dumpfile("waveform.vcd");
        $dumpvars(0, my_waveform_gem_tb);

        reset = 1;
        control = 9'sd0;
        #10;
        reset = 0;

        // Start test 1 - 3.29 MHz
        control = 9'sd8;
        #800;   // Roughly 2 time periods

        // Start test 2 - -1.35 MHz
        control = -9'sd3;
        #2000;  // Roughly 2 time periods

        // Start test 3 - 0.43 MHz
        control = 9'sd1;
        #5000; // Roughly 2 time periods

        $finish;

    end

    // Store the sine and cosine outputs to the file
    always @(posedge clk) begin
        if (!reset) begin
            $fwrite(outfile, "%d, %d\n", my_sine_out, my_cosine_out);
        end
    end



    endmodule








