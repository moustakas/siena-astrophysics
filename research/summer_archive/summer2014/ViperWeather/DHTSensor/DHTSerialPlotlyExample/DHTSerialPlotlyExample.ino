#include <SPI.h>
#include "DHT.h"


// DHT Sensor Setup
#define DHTPIN 2 // We have connected the DHT to Digital Pin 2
#define DHTTYPE DHT22 // This is the type of DHT Sensor (Change it to DHT11 if you're using that model)
DHT dht(DHTPIN, DHTTYPE); // Initialize DHT object


void setup() {
  // Open serial communications and wait for port to open:
  Serial.begin(9600);
  dht.begin();
}

float h, t;

void loop() {
      h = dht.readHumidity();
      t = dht.readTemperature();
      Serial.print(h);
      Serial.print("\t");
      Serial.println(t);
      delay(10000);
}
