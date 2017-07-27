using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using System.Collections.Generic;

public class GeodesicWorld : MonoBehaviour {

	public Vertex vertPrefab;
	//public Surface surface; // use for viewing in editor
	public List<Vertex> tiles;
	public int hexA, hexB; // defines the geodesic type
	public float radius;
	public Text inputDimX, inputDimY;

	private float _edgeErrorValue = .001f;
	// Use this for initialization
	private void Start() {
		CreateWorld();
	}
	
	public void CreateWorld() {
		for (int i = tiles.Count - 1; i >= 0; i--) {
			GameObject.DestroyImmediate(tiles[i].gameObject);
		}

		if (inputDimX.text != "" && inputDimY.text != "") {
			hexA = int.Parse(inputDimX.text);
			hexB = int.Parse(inputDimY.text);
		}
		tiles = GeodesicGenerator.GenerateGeodesic(hexA, hexB, vertPrefab);
		foreach (Vertex t in tiles) {
			t.transform.SetParent(transform);
		}
	}

	// Update is called once per frame
	private void Update() {

	}

	private void OnGUI() {
		
	}
}
