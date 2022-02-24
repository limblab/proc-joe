/***************************************************************************
  This is a library for the LPS22HB Absolute Digital Barometer

  Designed to work with all kinds of LPS22HB Breakout Boards

  These sensors use I2C, 2 pins are required to interface, as this :
	VDD to 3.3V DC
	SCL to A5
	SDA to A4
	GND to common groud 

  Written by Adrien Chapelet for Iotech
 ***************************************************************************/

#include <Wire.h>
#include <IO_LPS22HB.h>


IO_LPS22HB lps22hb;

void setup()
{

  pinMode(A4,INPUT);
  pinMode(A5,INPUT);
  
	Serial.begin(9600);
	
	lps22hb.begin(0x5C);

}

void loop()
{
	Serial.println(lps22hb.readPressure());
	delay(100);
}

