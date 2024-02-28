import * as THREE from "three";
import * as PHY from "simplePhysics";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

let renderer, scene, camera;
let world = {
  x: 80,
  z: 80,
};
let agentData = [];
const RADIUS = 1;
const defaultAgentMaterial = new THREE.MeshLambertMaterial({
  color: 0x00ff00,
});

init();
render();

function getVelocity() {
  return Math.random() - 0.5;
}

function init() {
  // renderer
  renderer = new THREE.WebGLRenderer();
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap; //
  renderer.setSize(window.innerWidth, window.innerHeight);
  document.body.appendChild(renderer.domElement);
  // scene
  scene = new THREE.Scene();
  // camera
  camera = new THREE.PerspectiveCamera(
    45,
    window.innerWidth / window.innerHeight,
    1,
    1000
  );

  camera.position.set(-67.26, 54.07, -3.77);
  camera.rotation.order = "YXZ";
  camera.rotation.y = -1.6267;
  camera.rotation.x = -0.46;

  // controls
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.addEventListener("change", render);
  controls.enableZoom = false;
  controls.enablePan = false;
  controls.maxPolarAngle = Math.PI / 2;

  // light
  const light = new THREE.PointLight(0xffffff, 0.9, 0, 100000);
  light.position.set(0, 50, 120);
  light.castShadow = true;
  light.shadow.mapSize.width = 512; // default
  light.shadow.mapSize.height = 512; // default
  light.shadow.camera.near = 0.5; // default
  light.shadow.camera.far = 5000; // default

  const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
  directionalLight.castShadow = true;
  directionalLight.position.set(-5, 20, 4);
  directionalLight.target.position.set(9, 0, -9);
  directionalLight.shadow.camera.left *= 9;
  directionalLight.shadow.camera.right *= 9;
  directionalLight.shadow.camera.top *= 9;
  directionalLight.shadow.camera.bottom *= 9;
  scene.add(directionalLight);

  // axes
  scene.add(new THREE.AxesHelper(40));
  const loader = new THREE.TextureLoader();
  const texture = loader.load("resources/OIP.jpg");
  texture.wrapS = THREE.RepeatWrapping;
  texture.wrapT = THREE.RepeatWrapping;
  texture.magFilter = THREE.NearestFilter;
  const repeats = 40 / 32;
  texture.repeat.set(repeats, repeats);

  // grid
  const geometry = new THREE.PlaneGeometry(world.x, world.z, 10, 10);
  const material = new THREE.MeshPhongMaterial({
    map: texture,
    //side: THREE.DoubleSide,
  });
  const grid = new THREE.Mesh(geometry, material);
  grid.castShadow = true; //default is false
  grid.receiveShadow = true; //default
  grid.rotation.order = "YXZ";
  grid.rotation.y = -Math.PI / 2;
  grid.rotation.x = -Math.PI / 2;
  scene.add(grid);

  let i = 0;
  let deltaSpacing = 5;
  let startX, startY, goalX, goalY;
  startX = -25;
  goalX = -25;
  startY = -30;
  goalY = 30;
  while (i < 10) {
    agentData.push({
      index: i,
      x: startX + i * deltaSpacing,
      y: 2.0,
      z: startY,
      goal_x: goalX + i * deltaSpacing,
      goal_y: 0.0,
      goal_z: goalY,
      vx: getVelocity(),
      vy: 0.0,
      vz: getVelocity(),
      radius: RADIUS,
      colliding: false,
    });
    agentData.push({
      index: i + 1,
      x: goalX + i * deltaSpacing + RADIUS,
      y: 2.0,
      z: goalY,
      goal_x: startX + i * deltaSpacing + RADIUS,
      goal_y: 0.0,
      goal_z: startY,
      vx: getVelocity(),
      vy: 0.0,
      vz: getVelocity(),
      radius: RADIUS,
      colliding: false,
    });
    i += 1;
  }

  let agnetGeometry, agentMaterial, agent;

  agentData.forEach(function (item) {
    agnetGeometry = new THREE.CylinderGeometry(item.radius, 1, 4, 16);
    agentMaterial = new THREE.MeshLambertMaterial({
      color: 0x00ff00,
    });
    agent = new THREE.Mesh(agnetGeometry, agentMaterial);
    agent.castShadow = true;
    agent.receiveShadow = true;
    scene.add(agent);
    item.agent = agent;
  });
}

function render() {
  renderer.render(scene, camera);
}

function animate() {
  requestAnimationFrame(animate);
  PHY.step(RADIUS, agentData, world);
  agentData.forEach(function (member) {
    member.agent.position.x = member.x;
    member.agent.position.y = member.y;
    member.agent.position.z = member.z;
    member.agent.material = defaultAgentMaterial;
    /*  -------------------------------------  */
    /*  Optional TODO for debugging collisions */
    /*  -------------------------------------  */
  });
  renderer.render(scene, camera);
}

animate();
