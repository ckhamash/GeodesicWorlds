using UnityEngine;
using System.Collections;

[RequireComponent (typeof (Camera))]
public class GlobeCam : MonoBehaviour {

	GeodesicWorld globe;
	float distance = 10;
	private Vector3 guideAngles; //

	// Use this for initialization
	void Start () {
		globe = transform.parent.GetComponent<GeodesicWorld>();
		transform.position = new Vector3(distance, 0, 0);
		transform.LookAt(globe.transform, Vector3.up);
	}
	
	// Update is called once per frame
	void Update () {

		distance -= Input.mouseScrollDelta.y;
		Vector3 rotateAngles = new Vector3(0, 0, 0);
		
		if (Input.GetKey(KeyCode.W)) {
			rotateAngles += new Vector3(5, 0, 0);
		}
		if (Input.GetKey(KeyCode.S)) {
			rotateAngles -= new Vector3(5, 0, 0);
		}
		if (Input.GetKey(KeyCode.A)) {
			rotateAngles += new Vector3(0, 5, 0);
		}
		if (Input.GetKey(KeyCode.D)) {
			rotateAngles -= new Vector3(0, 5, 0);
		}
		if (Input.GetKey(KeyCode.Q)) {
			rotateAngles += new Vector3(0, 0, 5);
		}
		if (Input.GetKey(KeyCode.E)) {
			rotateAngles -= new Vector3(0, 0, 5);
		}
		// use Quaternion.AngleAxis(angle, axis)
		transform.position =
			Quaternion.AngleAxis(rotateAngles.x, transform.right) *
			Quaternion.AngleAxis(rotateAngles.y, transform.up) *
			transform.position;

		transform.rotation = Quaternion.AngleAxis(rotateAngles.z, transform.forward) * transform.rotation;

		// scrolling/zoom
		for (int i = 0; i < Mathf.Abs(Input.mouseScrollDelta.y); i++) {
			if (Mathf.Sign(Input.mouseScrollDelta.y) > 0)
				transform.position *= 1.05f;
			else
				transform.position /= 1.05f;
		}

		transform.LookAt(globe.transform, transform.up);
	}

	//
	public void Initialize() {
		transform.position = new Vector3(globe.radius + distance, 0, 0);
		transform.LookAt(globe.transform, Vector3.up);
	}
}
