#include "Wire.h"

#define I2CAddress 0x27
long aa, bb, cc, dd, ee;
float temperature = 0.00;
float pressure = 0.00;


void setup() {
  // put your setup code here, to run once:
  pinMode(13, OUTPUT);
  digitalWrite(13, HIGH);        //switches sensor on
  Serial.begin(9600);  // start serial for output
  Wire.begin();
}


void NPA201ReadData()
{
  // Initiate Comms to device, initiate measure and read 5 bytes of data
  Wire.beginTransmission(I2CAddress);
  Wire.write(0xAC);  Wire.write(0);
  Wire.endTransmission();
  delay(21);
  Wire.requestFrom(I2CAddress, 5);
  aa = Wire.read();
  bb = Wire.read();
  cc = Wire.read();
  dd = Wire.read();
  ee = Wire.read();

  //Pressure Value = BridgeData/65535 * (1260-260) + 260 (hPa)
  //Temperature_Value = TempData/65535 * (85+40) - 40   (°C)
  // additional calculations to make values reasonable based on accuracy

  pressure =  (float) ((bb << 8) | cc) / 65535.0 * 1000.0 + 260;
  //pressure = round(pressure);

  temperature =  (float) ((dd << 8) | ee) / 65535 * 125 - 40;
  temperature = round(temperature * 10);
  temperature = temperature / 10;
}

void loop() {
  // put your main code here, to run repeatedly:
  NPA201ReadData();
  
  Serial.println((pressure));

  delay(10);
}
