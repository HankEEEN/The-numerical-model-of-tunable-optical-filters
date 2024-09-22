module Cosine_LUT (
    input [7:0] phase,
    output reg signed [7:0] cosine_value
);

    always @(*) begin 
        case (phase)
            8'd0: cosine_value = 8'b01111111;
            8'd1: cosine_value = 8'b01111111;
            8'd2: cosine_value = 8'b01111111;
            8'd3: cosine_value = 8'b01111111;
            8'd4: cosine_value = 8'b01111111;
            8'd5: cosine_value = 8'b01111111;
            8'd6: cosine_value = 8'b01111111;
            8'd7: cosine_value = 8'b01111110;
            8'd8: cosine_value = 8'b01111110;
            8'd9: cosine_value = 8'b01111101;
            8'd10: cosine_value = 8'b01111100;
            8'd11: cosine_value = 8'b01111011;
            8'd12: cosine_value = 8'b01111010;
            8'd13: cosine_value = 8'b01111010;
            8'd14: cosine_value = 8'b01111001;
            8'd15: cosine_value = 8'b01110111;
            8'd16: cosine_value = 8'b01110110;
            8'd17: cosine_value = 8'b01110101;
            8'd18: cosine_value = 8'b01110100;
            8'd19: cosine_value = 8'b01110010;
            8'd20: cosine_value = 8'b01110001;
            8'd21: cosine_value = 8'b01101111;
            8'd22: cosine_value = 8'b01101110;
            8'd23: cosine_value = 8'b01101100;
            8'd24: cosine_value = 8'b01101010;
            8'd25: cosine_value = 8'b01101001;
            8'd26: cosine_value = 8'b01100111;
            8'd27: cosine_value = 8'b01100101;
            8'd28: cosine_value = 8'b01100011;
            8'd29: cosine_value = 8'b01100001;
            8'd30: cosine_value = 8'b01011111;
            8'd31: cosine_value = 8'b01011101;
            8'd32: cosine_value = 8'b01011011;
            8'd33: cosine_value = 8'b01011000;
            8'd34: cosine_value = 8'b01010110;
            8'd35: cosine_value = 8'b01010100;
            8'd36: cosine_value = 8'b01010001;
            8'd37: cosine_value = 8'b01001111;
            8'd38: cosine_value = 8'b01001100;
            8'd39: cosine_value = 8'b01001010;
            8'd40: cosine_value = 8'b01000111;
            8'd41: cosine_value = 8'b01000100;
            8'd42: cosine_value = 8'b01000010;
            8'd43: cosine_value = 8'b00111111;
            8'd44: cosine_value = 8'b00111100;
            8'd45: cosine_value = 8'b00111010;
            8'd46: cosine_value = 8'b00110111;
            8'd47: cosine_value = 8'b00110100;
            8'd48: cosine_value = 8'b00110001;
            8'd49: cosine_value = 8'b00101110;
            8'd50: cosine_value = 8'b00101011;
            8'd51: cosine_value = 8'b00101000;
            8'd52: cosine_value = 8'b00100101;
            8'd53: cosine_value = 8'b00100010;
            8'd54: cosine_value = 8'b00011111;
            8'd55: cosine_value = 8'b00011100;
            8'd56: cosine_value = 8'b00011001;
            8'd57: cosine_value = 8'b00010110;
            8'd58: cosine_value = 8'b00010011;
            8'd59: cosine_value = 8'b00010000;
            8'd60: cosine_value = 8'b00001101;
            8'd61: cosine_value = 8'b00001001;
            8'd62: cosine_value = 8'b00000110;
            8'd63: cosine_value = 8'b00000011;
            8'd64: cosine_value = 8'b00000000;
            8'd65: cosine_value = 8'b11111101;
            8'd66: cosine_value = 8'b11111010;
            8'd67: cosine_value = 8'b11110111;
            8'd68: cosine_value = 8'b11110011;
            8'd69: cosine_value = 8'b11110000;
            8'd70: cosine_value = 8'b11101101;
            8'd71: cosine_value = 8'b11101010;
            8'd72: cosine_value = 8'b11100111;
            8'd73: cosine_value = 8'b11100100;
            8'd74: cosine_value = 8'b11100001;
            8'd75: cosine_value = 8'b11011110;
            8'd76: cosine_value = 8'b11011011;
            8'd77: cosine_value = 8'b11011000;
            8'd78: cosine_value = 8'b11010101;
            8'd79: cosine_value = 8'b11010010;
            8'd80: cosine_value = 8'b11001111;
            8'd81: cosine_value = 8'b11001100;
            8'd82: cosine_value = 8'b11001001;
            8'd83: cosine_value = 8'b11000110;
            8'd84: cosine_value = 8'b11000100;
            8'd85: cosine_value = 8'b11000001;
            8'd86: cosine_value = 8'b10111110;
            8'd87: cosine_value = 8'b10111100;
            8'd88: cosine_value = 8'b10111001;
            8'd89: cosine_value = 8'b10110110;
            8'd90: cosine_value = 8'b10110100;
            8'd91: cosine_value = 8'b10110001;
            8'd92: cosine_value = 8'b10101111;
            8'd93: cosine_value = 8'b10101100;
            8'd94: cosine_value = 8'b10101010;
            8'd95: cosine_value = 8'b10101000;
            8'd96: cosine_value = 8'b10100101;
            8'd97: cosine_value = 8'b10100011;
            8'd98: cosine_value = 8'b10100001;
            8'd99: cosine_value = 8'b10011111;
            8'd100: cosine_value = 8'b10011101;
            8'd101: cosine_value = 8'b10011011;
            8'd102: cosine_value = 8'b10011001;
            8'd103: cosine_value = 8'b10010111;
            8'd104: cosine_value = 8'b10010110;
            8'd105: cosine_value = 8'b10010100;
            8'd106: cosine_value = 8'b10010010;
            8'd107: cosine_value = 8'b10010001;
            8'd108: cosine_value = 8'b10001111;
            8'd109: cosine_value = 8'b10001110;
            8'd110: cosine_value = 8'b10001100;
            8'd111: cosine_value = 8'b10001011;
            8'd112: cosine_value = 8'b10001010;
            8'd113: cosine_value = 8'b10001001;
            8'd114: cosine_value = 8'b10000111;
            8'd115: cosine_value = 8'b10000110;
            8'd116: cosine_value = 8'b10000110;
            8'd117: cosine_value = 8'b10000101;
            8'd118: cosine_value = 8'b10000100;
            8'd119: cosine_value = 8'b10000011;
            8'd120: cosine_value = 8'b10000010;
            8'd121: cosine_value = 8'b10000010;
            8'd122: cosine_value = 8'b10000001;
            8'd123: cosine_value = 8'b10000001;
            8'd124: cosine_value = 8'b10000001;
            8'd125: cosine_value = 8'b10000000;
            8'd126: cosine_value = 8'b10000000;
            8'd127: cosine_value = 8'b10000000;
            8'd128: cosine_value = 8'b10000000;
            8'd129: cosine_value = 8'b10000000;
            8'd130: cosine_value = 8'b10000000;
            8'd131: cosine_value = 8'b10000000;
            8'd132: cosine_value = 8'b10000001;
            8'd133: cosine_value = 8'b10000001;
            8'd134: cosine_value = 8'b10000001;
            8'd135: cosine_value = 8'b10000010;
            8'd136: cosine_value = 8'b10000010;
            8'd137: cosine_value = 8'b10000011;
            8'd138: cosine_value = 8'b10000100;
            8'd139: cosine_value = 8'b10000101;
            8'd140: cosine_value = 8'b10000110;
            8'd141: cosine_value = 8'b10000110;
            8'd142: cosine_value = 8'b10000111;
            8'd143: cosine_value = 8'b10001001;
            8'd144: cosine_value = 8'b10001010;
            8'd145: cosine_value = 8'b10001011;
            8'd146: cosine_value = 8'b10001100;
            8'd147: cosine_value = 8'b10001110;
            8'd148: cosine_value = 8'b10001111;
            8'd149: cosine_value = 8'b10010001;
            8'd150: cosine_value = 8'b10010010;
            8'd151: cosine_value = 8'b10010100;
            8'd152: cosine_value = 8'b10010110;
            8'd153: cosine_value = 8'b10010111;
            8'd154: cosine_value = 8'b10011001;
            8'd155: cosine_value = 8'b10011011;
            8'd156: cosine_value = 8'b10011101;
            8'd157: cosine_value = 8'b10011111;
            8'd158: cosine_value = 8'b10100001;
            8'd159: cosine_value = 8'b10100011;
            8'd160: cosine_value = 8'b10100101;
            8'd161: cosine_value = 8'b10101000;
            8'd162: cosine_value = 8'b10101010;
            8'd163: cosine_value = 8'b10101100;
            8'd164: cosine_value = 8'b10101111;
            8'd165: cosine_value = 8'b10110001;
            8'd166: cosine_value = 8'b10110100;
            8'd167: cosine_value = 8'b10110110;
            8'd168: cosine_value = 8'b10111001;
            8'd169: cosine_value = 8'b10111100;
            8'd170: cosine_value = 8'b10111110;
            8'd171: cosine_value = 8'b11000001;
            8'd172: cosine_value = 8'b11000100;
            8'd173: cosine_value = 8'b11000110;
            8'd174: cosine_value = 8'b11001001;
            8'd175: cosine_value = 8'b11001100;
            8'd176: cosine_value = 8'b11001111;
            8'd177: cosine_value = 8'b11010010;
            8'd178: cosine_value = 8'b11010101;
            8'd179: cosine_value = 8'b11011000;
            8'd180: cosine_value = 8'b11011011;
            8'd181: cosine_value = 8'b11011110;
            8'd182: cosine_value = 8'b11100001;
            8'd183: cosine_value = 8'b11100100;
            8'd184: cosine_value = 8'b11100111;
            8'd185: cosine_value = 8'b11101010;
            8'd186: cosine_value = 8'b11101101;
            8'd187: cosine_value = 8'b11110000;
            8'd188: cosine_value = 8'b11110011;
            8'd189: cosine_value = 8'b11110111;
            8'd190: cosine_value = 8'b11111010;
            8'd191: cosine_value = 8'b11111101;
            8'd192: cosine_value = 8'b00000000;
            8'd193: cosine_value = 8'b00000011;
            8'd194: cosine_value = 8'b00000110;
            8'd195: cosine_value = 8'b00001001;
            8'd196: cosine_value = 8'b00001101;
            8'd197: cosine_value = 8'b00010000;
            8'd198: cosine_value = 8'b00010011;
            8'd199: cosine_value = 8'b00010110;
            8'd200: cosine_value = 8'b00011001;
            8'd201: cosine_value = 8'b00011100;
            8'd202: cosine_value = 8'b00011111;
            8'd203: cosine_value = 8'b00100010;
            8'd204: cosine_value = 8'b00100101;
            8'd205: cosine_value = 8'b00101000;
            8'd206: cosine_value = 8'b00101011;
            8'd207: cosine_value = 8'b00101110;
            8'd208: cosine_value = 8'b00110001;
            8'd209: cosine_value = 8'b00110100;
            8'd210: cosine_value = 8'b00110111;
            8'd211: cosine_value = 8'b00111010;
            8'd212: cosine_value = 8'b00111100;
            8'd213: cosine_value = 8'b00111111;
            8'd214: cosine_value = 8'b01000010;
            8'd215: cosine_value = 8'b01000100;
            8'd216: cosine_value = 8'b01000111;
            8'd217: cosine_value = 8'b01001010;
            8'd218: cosine_value = 8'b01001100;
            8'd219: cosine_value = 8'b01001111;
            8'd220: cosine_value = 8'b01010001;
            8'd221: cosine_value = 8'b01010100;
            8'd222: cosine_value = 8'b01010110;
            8'd223: cosine_value = 8'b01011000;
            8'd224: cosine_value = 8'b01011011;
            8'd225: cosine_value = 8'b01011101;
            8'd226: cosine_value = 8'b01011111;
            8'd227: cosine_value = 8'b01100001;
            8'd228: cosine_value = 8'b01100011;
            8'd229: cosine_value = 8'b01100101;
            8'd230: cosine_value = 8'b01100111;
            8'd231: cosine_value = 8'b01101001;
            8'd232: cosine_value = 8'b01101010;
            8'd233: cosine_value = 8'b01101100;
            8'd234: cosine_value = 8'b01101110;
            8'd235: cosine_value = 8'b01101111;
            8'd236: cosine_value = 8'b01110001;
            8'd237: cosine_value = 8'b01110010;
            8'd238: cosine_value = 8'b01110100;
            8'd239: cosine_value = 8'b01110101;
            8'd240: cosine_value = 8'b01110110;
            8'd241: cosine_value = 8'b01110111;
            8'd242: cosine_value = 8'b01111001;
            8'd243: cosine_value = 8'b01111010;
            8'd244: cosine_value = 8'b01111010;
            8'd245: cosine_value = 8'b01111011;
            8'd246: cosine_value = 8'b01111100;
            8'd247: cosine_value = 8'b01111101;
            8'd248: cosine_value = 8'b01111110;
            8'd249: cosine_value = 8'b01111110;
            8'd250: cosine_value = 8'b01111111;
            8'd251: cosine_value = 8'b01111111;
            8'd252: cosine_value = 8'b01111111;
            8'd253: cosine_value = 8'b01111111;
            8'd254: cosine_value = 8'b01111111;
            8'd255: cosine_value = 8'b01111111;
            default: cosine_value = 8'd0;
        endcase
    end
endmodule



module my_waveform_gem (
    input clk,
    input reset,
    input signed [8:0] control,
    output [7:0] my_sine_out,
    output [7:0] my_cosine_out
);

    reg signed [8:0] phase_change;
    

    always @(posedge clk) begin
        if (reset) 
            phase_change <= 9'sd0;
        else
            phase_change <= phase_change + control;
    end

    // Instantiation for cosine LUT output
    wire [7:0] lut_entry = phase_change [7:0];

    Cosine_LUT cos_ut (
        .phase(lut_entry),
        .cosine_value(my_cosine_out)
    );

    // Instantiation for sine LUT output
    Cosine_LUT sin_ut (
        .phase(lut_entry + 8'd64),  // 90-degree phase shift equivalent in terms of entry number
        .cosine_value(my_sine_out)
    );


endmodule


module my_controller (
    input clk, 
    input reset,
    input [7:0] period,
    output reg signed [8:0] control
);


    reg [13:0] counter;
    reg [13:0] half_period; // The number of the clock cycles
    reg signed [8:0] control_hf;
    reg signed [8:0] control_lf;
    reg switch; // Help to make switch between the two frequency


    // Initialize the control signal for each frequency (high for hf and low for lf)
    initial begin
        control_hf = 9'sd6;
        control_lf = 9'sd1;

    end


    // Alternates the frequency between 2.3 MHz and 0.3 MHz with the input period
    // reset cases - reset and period change
    always @(posedge clk) begin
        if (reset) begin
            counter <= 14'd0;
            switch <= 1'b0;
            control <= control_hf;
            half_period <= period * 8'd50;
        end else begin
        // Alternate the frequency when counter reaches the number of cycles for each half period
        if (counter >= half_period - 1) begin
            counter <= 14'd0;
            switch <= ~switch;
            control <= (~switch) ? control_lf : control_hf;
        end else begin
            counter <= counter + 14'd1;
        end
        // $display("counter number is:\t", counter);
        // $display("switch is:\t", switch);
        end
    end

endmodule


