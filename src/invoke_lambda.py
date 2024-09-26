from marshal import loads

import boto3
import json
from dotenv import load_dotenv

load_dotenv()

session = boto3.Session()
lambda_client = session.client('lambda')

event = {
    "names": ["Alice", "Bob", "Charlie"]
}

response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = json.loads(response['Payload'].read())

print(response_payload)