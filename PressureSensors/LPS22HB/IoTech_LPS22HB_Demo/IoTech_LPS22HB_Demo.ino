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
	Serial.println("IoTech LPS22HB Arduino Test");
	
	lps22hb.begin(0x5C);

	Serial.print("Who Am I? 0x");
	Serial.print(lps22hb.whoAmI(), HEX);
	Serial.println(" (expected: 0xB1)");
}

void loop()
{
	Serial.print("P=");
	Serial.print(lps22hb.readPressure());
	Serial.print(" mbar, T=");
	Serial.print((lps22hb.readTemperature()-42.5f)*480.0f);
	Serial.println("C");
	delay(300);
}

