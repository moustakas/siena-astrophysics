#include <SPI.h>
#include <Ethernet.h>
#include "plotly_streaming_ethernet.h"
#include "DHT.h"

// Sign up to plotly here: https://plot.ly
// View your API key and streamtokens here: https://plot.ly/settings
#define nTraces 2
// View your tokens here: https://plot.ly/settings
// Supply as many tokens as data traces
// e.g. if you want to ploty A0 and A1 vs time, supply two tokens
char *tokens[nTraces] = {"kijrfl3wbp", "m7ngmur9bl"};
// arguments: username, api key, streaming token, filename
plotly graph("janedoe", "pword", tokens, "DHT", nTraces);

// DHT Sensor Setup
#define DHTPIN 2 // We have connected the DHT to Digital Pin 2
#define DHTTYPE DHT22 // This is the type of DHT Sensor (Change it to DHT11 if you're using that model)
DHT dht(DHTPIN, DHTTYPE); // Initialize DHT object

byte mac[] = {0x00, 0xAA, 0xBB, 0xCC, 0xDE, 0x02  };//the mac address of the ethernet shield (made it up!)
byte my_ip[] = {10,0,1,231 }; // the ip of your ethernet shield

void startEthernet(){
    Serial.println("... Initializing ethernet");
    if(Ethernet.begin(mac) == 0){
        Serial.println("... Failed to configure Ethernet using DHCP");
        // no point in carrying on, so do nothing forevermore:
        // try to congifure using IP address instead of DHCP:
        Ethernet.begin(mac, my_ip);
    }
    Serial.println("... Done initializing ethernet");
    delay(1000);
}


void setup() {
  graph.maxpoints = 100;
  // Open serial communications and wait for port to open:
  Serial.begin(9600);
  while (!Serial) {
    ; // wait for serial port to connect. Needed for Leonardo only
  }

  startEthernet();

  bool success;
  success = graph.init();
  if(!success){while(true){}}
  graph.openStream();
  dht.begin();
}

float h, t;

void loop() {
      h = dht.readHumidity();
      t = dht.readTemperature();
      graph.plot(millis(), t, tokens[0]);
      delay(1000);
      graph.plot(millis(), h, tokens[1]);
      delay(10000);
}
