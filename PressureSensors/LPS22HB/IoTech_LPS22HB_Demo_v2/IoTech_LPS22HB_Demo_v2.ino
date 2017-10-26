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


/* 
 *  
 *  DEFINE SENSOR CLASS
 *  
 */
 
#define _IO_LPS22HB_h
#include <Arduino.h>

#define LPS22HB_WHO_AM_I  0X0F //Who am I
#define LPS22HB_RES_CONF  0X1A //Resolution
#define LPS22HB_CTRL_REG1 0X10
#define LPS22HB_CTRL_REG2 0X11
#define LPS22HB_STATUS_REG  0X27
#define LPS22HB_PRES_OUT_XL 0X28 //LSB
#define LPS22HB_PRES_OUT_L  0X29
#define LPS22HB_PRES_OUT_H  0X2A //MSB
#define LPS22HB_TEMP_OUT_L  0X2B //LSB
#define LPS22HB_TEMP_OUT_H  0X2C //MSB



class IO_LPS22HB {
  public:
    IO_LPS22HB();
  
    void begin(uint8_t address);
  
    uint8_t whoAmI();
    float readTemperature();
  
    float readPressure();
    uint32_t readPressureUI();
    uint32_t readPressureRAW();
  
  private:
    uint8_t _address;
    uint8_t read(uint8_t reg);
    void write(uint8_t reg, uint8_t data);
    uint8_t status(uint8_t data);
};

/* 
 *  
 *  CLASS FUNCTIONS BELOW
 *  
 */
 
IO_LPS22HB::IO_LPS22HB() // constructor
{
}

void IO_LPS22HB::begin(uint8_t address) {
  _address = address;
  Wire.begin();
  write(LPS22HB_RES_CONF, 0x0); // resolution: temp=32, pressure=128
  write(LPS22HB_CTRL_REG1, 0b01010000); // 75Hz continuous streaming
}

byte IO_LPS22HB::whoAmI() {
  Wire.beginTransmission(_address);
  Wire.write(LPS22HB_WHO_AM_I);
  Wire.endTransmission();
  Wire.requestFrom(_address, 1);
  return Wire.read();
}

float IO_LPS22HB::readPressure() {
  write(LPS22HB_CTRL_REG2, 0x1);

  if (status(0x1) < 0)
    return 1.23;
  //delay(50);
  uint8_t pressOutH = read(LPS22HB_PRES_OUT_H);
  uint8_t pressOutL = read(LPS22HB_PRES_OUT_L);
  uint8_t pressOutXL = read(LPS22HB_PRES_OUT_XL);

  long val = ( ((long)pressOutH << 24) | ((long)pressOutL << 16) | ((long)pressOutXL << 8)) >> 8;
  //if (val == 1.00) readPressure();
  return val/4096.0f;
}

uint32_t IO_LPS22HB::readPressureRAW() {
  write(LPS22HB_CTRL_REG2, 0x1);

  if (status(0x1) < 0)
    return 123;
  //delay(50);
  uint8_t pressOutH = read(LPS22HB_PRES_OUT_H);
  uint8_t pressOutL = read(LPS22HB_PRES_OUT_L);
  uint8_t pressOutXL = read(LPS22HB_PRES_OUT_XL);

  int32_t val = ((pressOutH << 24) | (pressOutL << 16) | (pressOutXL << 8));
  val >> 8;
  val=val+0x400000;
  //if (val == 1.00) readPressure();
  return (uint32_t)val;
}

uint32_t IO_LPS22HB::readPressureUI() {
  write(LPS22HB_CTRL_REG2, 0x1);

  if (status(0x1) < 0)
    return 1.23;
  //delay(50);
  uint8_t pressOutH = read(LPS22HB_PRES_OUT_H);
  uint8_t pressOutL = read(LPS22HB_PRES_OUT_L);
  uint8_t pressOutXL = read(LPS22HB_PRES_OUT_XL);

  uint32_t val = ((pressOutH << 24) | (pressOutL << 16) | (pressOutXL << 8)) >> 8;
  //if (val == 1.00) readPressure();
  return val/4096;
}

float IO_LPS22HB::readTemperature() {
  write(LPS22HB_CTRL_REG2, 0x1);
  if (status(0x2) < 0)
    return 4.56;

  uint8_t tempOutH = read(LPS22HB_TEMP_OUT_H);
  uint8_t tempOutL = read(LPS22HB_TEMP_OUT_L);

  int16_t val = tempOutH << 8 | tempOutL & 0xff;
  return 42.5f+val/480.0f;
}


uint8_t IO_LPS22HB::status(uint8_t status) {
  int count = 1000;
  uint8_t data = 0xff;
  do {
    data = read(LPS22HB_STATUS_REG);
    --count;
    if (count < 0)
      break;
  } while ((data & status) == 0);

  if (count < 0)
    return -1;
  else
    return 0;
}

uint8_t IO_LPS22HB::read(uint8_t reg) {
  Wire.beginTransmission(_address);
  Wire.write(reg);
  Wire.endTransmission();
  Wire.requestFrom(_address, 1);
  return Wire.read();
}

void IO_LPS22HB::write(uint8_t reg, uint8_t data) {
  Wire.beginTransmission(_address);
  Wire.write(reg);
  Wire.write(data);
  Wire.endTransmission();
}

// end define class functions


/* 


SETUP AND LOOP CODE BELOW 


*/


// create class object
IO_LPS22HB lps22hb = IO_LPS22HB();

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
}

