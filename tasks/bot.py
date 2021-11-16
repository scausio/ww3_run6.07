import json
import requests
from datetime import datetime
'''import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--message', '-m', help='message to send')

args = parser.parse_args()

TOKEN = '625361077:AAHytmtu-kIiSKXwwk5zuTZ4YdVy3pIKDM4'
URL = "https://api.telegram.org/bot%s/" % (TOKEN)
'''
def get_url(url):
    response = requests.get(url)
    content = response.content.decode("utf8")
    return content

def get_json_from_url(url):
    content = get_url(url)
    js = json.loads(content)
    return js

def get_updates(URL):
    url = URL + "getUpdates"
    js = get_json_from_url(url)
    return js

def get_last_chat_id_and_text(updates):
    num_updates = len(updates["result"])
    last_update = num_updates - 1
    text = updates["result"][last_update]["message"]["text"]
    chat_id = updates["result"][last_update]["message"]["chat"]["id"]
    return (text, chat_id)

def send_message(text):
    token = "625361077:AAHytmtu-kIiSKXwwk5zuTZ4YdVy3pIKDM4"
    address = "https://api.telegram.org/bot%s/" % (token)
    chat_id="106157864"
    url =  "{u}sendMessage?text={t}&chat_id={c}".format(u=address,t=text,c=chat_id)
    get_updates(address)
    get_url(url)


#send_message(args.message, '106157864')


