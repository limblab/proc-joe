float freq = 5;
int num_sweeps = 20;
int ttl_volt = 255;
int zero_volt = 0;

void setup() {
  // put your setup code here, to run once:
  pinMode(A0,OUTPUT);


// stim code
  for(int i = 0; i < num_sweeps; i++){
    // each sweep, sweep through the pulse widths
    for(int pw = 50; pw < 1000; pw = pw + 50){
      // set pin_out high
      analogWrite(A0,ttl_volt);
      // wait for pw microseconds
      delayMicroseconds(pw);
      // set pin_out low
      analogWrite(A0,zero_volt);
      // wait for 1/freq
      delay(1000/freq);
    }
    delay(100);
  }
  
}

void loop() {
 // do nothing
}
