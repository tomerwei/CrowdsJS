import * as THREE from "three";

export function distance(x1, y1, x2, y2) {
  return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// Function to calculate the angle between two vectors
function angleBetweenVectors_2(v1, v2) {
  // Calculate dot product
  const dotProduct = v1.x * v2.x + v1.z * v2.z;

  // Calculate magnitudes of vectors
  const magnitudeV1 = Math.sqrt(v1.x * v1.x + v1.z * v1.z);
  const magnitudeV2 = Math.sqrt(v2.x * v2.x + v2.z * v2.z);

  // Calculate angle in radians using dot product and magnitudes
  const angleRadians = Math.acos(dotProduct / ((magnitudeV1 * magnitudeV2)+0.0000001) );

  // Convert angle from radians to degrees
  const angleDegrees = angleRadians * 180 / Math.PI;

  return angleDegrees;
}


// Function to calculate the angle between two vectors
function calculateAngle(point1_x, point1_z, point2_x, point2_z, velocity_x, velocity_z){

  // Calculate vectors from point1 to point2 and from point1 to the velocity
  const vector1 = { x: point2_x - point1_x, z: point2_z - point1_z };
  const vector2 = { x: velocity_x - point1_x, z: velocity_z - point1_z };


  // Calculate the dot product of the two vectors
  const dotProduct = vector1.x * vector2.x + vector1.z * vector2.z;

  // Calculate the magnitudes of the vectors
  const magnitude1 = Math.sqrt(vector1.x * vector1.x + vector1.z * vector1.z);
  const magnitude2 = Math.sqrt(vector2.x * vector2.x + vector2.z * vector2.z);

  // Calculate the cosine of the angle
  const cosineTheta = dotProduct / (magnitude1 * magnitude2);

  // Calculate the angle in radians using the arccosine
  const angleRadians = Math.acos(cosineTheta);

  // Convert the angle to degrees
  const angleDegrees = angleRadians * (180 / Math.PI);

  return angleDegrees;
}

function degreeBetween(p1x, p1z, p2x, p2z){

  // Calculate the dot product of the two vectors
  const dotProduct = p1x * p2x + p1z * p2z;

  // Calculate the magnitude (length) of the vectors
  const magnitude1 = Math.sqrt(p1x * p1x + p1z * p1z);
  const magnitude2 = Math.sqrt(p2x * p2x + p2z * p2z);

  // Calculate the cosine of the angle between the vectors
  const cosTheta = dotProduct / (magnitude1 * magnitude2);

  // Calculate the angle in radians
  const angleRad = Math.acos(cosTheta);

  // Convert radians to degrees
  const angleDeg = angleRad * (180 / Math.PI);

  // Calculate the signed angle in degrees using the atan2 function
  const signedAngleDeg = Math.atan2(p2z, p2x) - Math.atan2(p1z, p2x);

  // Ensure the angle is between 0 and 360 degrees
  const normalizedAngleDeg = (signedAngleDeg >= 0) ? signedAngleDeg : (signedAngleDeg + 2 * Math.PI) * (180 / Math.PI);

return angleDeg;
}


// Function to get a vector direction deviated by a certain angle from a known normalized vector
function getNewVectorDeviatedFromAnother(knownVector, angle) {
  // Calculate the components of the known vector
  var x = knownVector.x;
  var y = knownVector.y;
  var z = knownVector.z;

  // Calculate the new vector direction using trigonometric functions
  var newX = x * Math.cos(angle) + (1 - Math.cos(angle)) * x + Math.sin(angle) * (y - z);
  var newY = y * Math.cos(angle) + (1 - Math.cos(angle)) * y + Math.sin(angle) * (z - x);
  var newZ = z * Math.cos(angle) + (1 - Math.cos(angle)) * z + Math.sin(angle) * (x - y);

  // Return the new vector direction as a normalized vector
  var magnitude = Math.sqrt(newX * newX + newY * newY + newZ * newZ);
  return { x: newX / magnitude, y: 0, z: newZ / magnitude };
}


function findNewPositionByProjecttion(predict_point, cur_point, direction_vector) {
  // get position vector from current and precidted position
  const position_vector = new THREE.Vector3().subVectors(predict_point, cur_point);

  // Project position_vector onto w_r
  const proj = position_vector.clone().projectOnVector(direction_vector); 

  // Calculate new predicted positionon on w_r vector direction.
  const p_new = cur_point.clone().add(proj);

  return p_new;  
}



export function step(RADIUS, sceneEntities, world, scene, customParams = {}) {
  const AGENTSIZE = RADIUS * 2;
  const epsilon = 0.0001;
  const timestep = 0.03;
  const ITERNUM = 1; // 3
  const agentLength = RADIUS;

  let C_TAU_MAX = 20;
  let C_TAO0 = 250; 

  const C_LONG_RANGE_STIFF = 0.04;  //previously used before experimenting for orientation
  const MAX_DELTA = 0.05;

  let angle_threshold = 0.1299;   //angle less or equal to 0.12 doesn't work.

  // collision functions
  function rotateLineSegment(x1, y1, x2, y2, r) {
    // Calculate the center of the line segment
    const centerX = (x1 + x2) / 2;
    const centerY = (y1 + y2) / 2;

    // Translate the line segment so that its center is at the origin
    const x1p = x1 - centerX;
    const y1p = y1 - centerY;
    const x2p = x2 - centerX;
    const y2p = y2 - centerY;

    // Rotate the line segment about the origin
    const cosR = Math.cos(r);
    const sinR = Math.sin(r);
    const x1r = x1p * cosR - y1p * sinR;
    const y1r = x1p * sinR + y1p * cosR;
    const x2r = x2p * cosR - y2p * sinR;
    const y2r = x2p * sinR + y2p * cosR;

    // Translate the line segment back to its original position
    const newX1 = x1r + centerX;
    const newY1 = y1r + centerY;
    const newX2 = x2r + centerX;
    const newY2 = y2r + centerY;

    // Return the new endpoints of the line segment
    return [newX1, newY1, newX2, newY2];
  }

  function ClosestPointOnLineSegment(A, B, Point) {
    const AB = B.clone().sub(A);
    const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
    return A.clone().add(
      AB.clone().multiplyScalar(Math.min(Math.max(t, 0), 1))
    );
  }

  function PointOnLineSegment(A, B, Point) {
    const AB = B.clone().sub(A);
    const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
    return A.clone().add(
        AB.clone().multiplyScalar(t)
    );
  }

  function is_colliding_torso(x11, y11, x12, y12, x21, y21, x22, y22) {
    // console.log(segments_distance(x11, y11, x12, y12, x21, y21, x22, y22));
    return (
      segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) < 2 * RADIUS
    );
  }

  function segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) {
    // distance between two segments in the plane:
    // one segment is (x11, y11) to (x12, y12)
    // the other is (x21, y21) to (x22, y22)

    if (segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22)) return 0;

    // try each of the 4 vertices w/the other segment
    let distances = [
      point_segment_distance(x11, y11, x21, y21, x22, y22),
      point_segment_distance(x12, y12, x21, y21, x22, y22),
      point_segment_distance(x21, y21, x11, y11, x12, y12),
      point_segment_distance(x22, y22, x11, y11, x12, y12),
    ];
    return Math.min(...distances);
  }

  function segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22) {
    // whether two segments in the plane intersect:
    // one segment is (x11, y11) to (x12, y12)
    // the other is (x21, y21) to (x22, y22)

    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;

    let delta = dx2 * dy1 - dy2 * dx1;
    if (delta === 0) return false; // parallel segments

    let s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta;
    let t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / -delta;

    return 0 <= s && s <= 1 && 0 <= t && t <= 1;
  }

  function point_segment_distance(px, py, x1, y1, x2, y2) {
    let dx = x2 - x1;
    let dy = y2 - y1;

    if (dx === 0 && dy === 0) {
      // the segment's just a point
      return Math.hypot(px - x1, py - y1);
    }

    // Calculate the t that minimizes the distance.
    let t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

    // See if this represents one of the segment's
    // end points or a point in the middle.
    if (t < 0) {
      dx = px - x1;
      dy = py - y1;
    } else if (t > 1) {
      dx = px - x2;
      dy = py - y2;
    } else {
      let near_x = x1 + t * dx;
      let near_y = y1 + t * dy;
      dx = px - near_x;
      dy = py - near_y;
    }
    return Math.hypot(dx, dy);
  }

  /*  -----------------------  */
  /*  TODO modify lines below  */
  /*  -----------------------  */
  function distanceConstraint(agent_i, agent_j, desiredDistance) {
    const agentCentroidDist = distance(
      agent_i.px,
      agent_i.pz,
      agent_j.px,
      agent_j.pz
    );
    const agentDist = agentCentroidDist - desiredDistance;
    const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
    const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
    const agent_i_scaler =
      ((0.1 * agent_i.invmass) / (agent_i.invmass + agent_j.invmass)) *
      agentDist;
    const agent_j_scaler =
      ((0.1 * agent_j.invmass) / (agent_i.invmass + agent_j.invmass)) *
      agentDist;
    if (Math.abs(agentDist) > epsilon) {
      agent_i.px += agent_i_scaler * dir_x;
      agent_i.pz += agent_i_scaler * dir_z;
      agent_j.px += -agent_j_scaler * dir_x;
      agent_j.pz += -agent_j_scaler * dir_z;
    }
  }

  function collisionConstraint(agent_i, agent_j) {
    const agentCentroidDist = distance(
      agent_i.px,
      agent_i.pz,
      agent_j.px,
      agent_j.pz
    );
    const agentDist = agentCentroidDist - AGENTSIZE;
    const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
    const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
    const agent_i_scaler =
      (agent_i.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
    const agent_j_scaler =
      (agent_j.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
    if (agentDist < 0) {
      agent_i.px += agent_i_scaler * dir_x;
      agent_i.pz += agent_i_scaler * dir_z;
      agent_j.px += -agent_j_scaler * dir_x;
      agent_j.pz += -agent_j_scaler * dir_z;
    }
  }


  function collisionConstraint_Capsule(best_i, best_j, p_best_i, p_best_j) {

    const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

    const agentDist = agentCentroidDist - AGENTSIZE;
    const dir_x = (p_best_j.x - p_best_i.x) / agentCentroidDist;
    const dir_z = (p_best_j.z - p_best_i.z) / agentCentroidDist;
    const agent_i_scaler = (1 / (1 + 1)) * agentDist;
    const agent_j_scaler = (1 / (1 + 1)) * agentDist;

      if (agentDist < 0) {

      sceneEntities[i].px += agent_i_scaler * dir_x;
      sceneEntities[i].pz += agent_i_scaler * dir_z;
      sceneEntities[j].px += -agent_j_scaler * dir_x;
      sceneEntities[j].pz += -agent_j_scaler * dir_z;
        
      sceneEntities[i].grad.dx += agent_i_scaler * dir_x;
      sceneEntities[i].grad.dz += agent_i_scaler * dir_z;
      sceneEntities[j].grad.dx += -agent_j_scaler * dir_x;
      sceneEntities[j].grad.dz += -agent_j_scaler * dir_z;
    }

    return [
      sceneEntities[i].grad.x,
      sceneEntities[i].grad.z,
        [dir_x],
        [dir_z]
    ];
  } 


  function agentVelocityPlanner() {
    sceneEntities.forEach(function (agent_i) {
      const distToGoal = distance(
        agent_i.x,
        agent_i.z,
        agent_i.goal_x,
        agent_i.goal_z
      );
      if (distToGoal > RADIUS) {
        const dir_x = (agent_i.goal_x - agent_i.x) / distToGoal;
        const dir_z = (agent_i.goal_z - agent_i.z) / distToGoal;
        agent_i.vx = agent_i.v_pref * dir_x;
        agent_i.vz = agent_i.v_pref * dir_z;
      }
      agent_i.vx = 0.9999 * agent_i.vx;
      agent_i.vz = 0.9999 * agent_i.vz;
    });
  }

  function clamp2D(vx,vy, maxValue) {
    const lengthV = Math.sqrt(vx * vx + vy * vy);
    if (lengthV > maxValue) {
      const mult = (maxValue / lengthV);
      vx *= mult;
      vy *= mult;
    }
    return {"x":vx, "y":vy}
  }

 
  //**************************************************** START ************************************************************ */
  // long-range collision avoidance for capsule shape.
  function longRangeConstraintCapsule(best_i, best_j,
                                      p_best_i, p_best_j,
                                      theta_i, theta_j,
                                      agent_i, agent_j,
                                      entity_i, entity_j,
                                      i = -1, j = -1) {

    const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

    const radius_init = 2 * AGENTSIZE;
    const radius_sq_init = radius_init * radius_init;
    let radius_sq = radius_sq_init;
    const dv_i = 1.;  // 1./delta_t;
    let delta_correction_i = {"x":0, "y":0};
    let delta_correction_j= {"x":0, "y":0};
    if (agentCentroidDist < radius_init) {
      radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
    }
    const v_x = (p_best_i.x - best_i.x) / timestep - (p_best_j.x - best_j.x) / timestep;
    const v_y = (p_best_i.z - best_i.z) / timestep - (p_best_j.z - best_j.z) / timestep;
    const x0 = best_i.x - best_j.x;
    const y0 = best_i.z - best_j.z;
    const v_sq = v_x * v_x + v_y * v_y;
    const x0_sq = x0 * x0;
    const y0_sq = y0 * y0;
    const x_sq = x0_sq + y0_sq;
    const a = v_sq;
    const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
    const b_sq = b * b;
    const c = x_sq - radius_sq;
    const d_sq = b_sq - a * c;
    const d = Math.sqrt(d_sq);
    const tao = (b - d) / a;
    // console.log("ttc in long range paper: " + tao);

    let lengthV;
    let grad_x_i;
    let grad_y_i;
    let grad_x_j;
    let grad_y_j;
    let s;

    if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {
      const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
      const tao_sq = c_tao * c_tao;

      grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
      grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
      grad_x_j = -grad_x_i;
      grad_y_j = -grad_y_i;
      
      // special case

      // rotation (facing angle) difference
      let facingDiff = Math.abs(theta_i - theta_j);
      //console.log( "theta_i: ", theta_i, ", theta_j: ", theta_j);

      // best points distance difference
      let projectedPoint_j = PointOnLineSegment(agent_j.real_base, agent_j.real_tip, best_i);
      let bestPointDiff = distance(projectedPoint_j.x, projectedPoint_j.z, best_j.x, best_j.z);

      // if facing direction on the same line AND the best points are exactly facing with each other
      // adding gradient value
      // if (bestPointDiff <= 0.01 &&  (facingDiff < 0.01 || facingDiff - Math.PI) <0.01 ){
      // if (bestPointDiff <= 2.0 &&  (facingDiff < 0.1 || facingDiff - Math.PI) <0.1 ){      
      if (bestPointDiff <= 4.0 &&  (facingDiff < 0.015 || facingDiff - Math.PI) <0.015 ){  //used this
        
        // grad_y_i = signNoP(grad_y_i) * Math.abs(grad_x_i);
        grad_y_i = signNoP(grad_y_i) * (Math.abs(grad_x_i)/1.2);   //used this
        // grad_y_i = signNoP(grad_y_i) * (  Math.abs(grad_x_i)/1.4  ); 
        // grad_y_i = signNoP(grad_y_i) * (Math.abs(grad_x_i) * 1.2);

        grad_y_j = -grad_y_i;
      }
    
      const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
      s = stiff * tao_sq / (0.5 * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + 0.5 * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed
      // console.log()

      delta_correction_i = clamp2D(s * 0.5 * grad_x_i,
          s * 0.5 * grad_y_i,
          MAX_DELTA);
      delta_correction_j = clamp2D(s * 0.5 * grad_x_j,
          s * 0.5 * grad_y_j,
          MAX_DELTA);
    }else {
      grad_x_i = 0;
      grad_y_i = 0;
      grad_x_j = 0;
      grad_y_j = 0;
      s=0;
    }

    // return tao;
    return [
        delta_correction_i,
        delta_correction_j,
        [grad_x_i, grad_y_i],
        [grad_x_j, grad_y_j],
        s
    ];
  }

  function getBestPointWithWall(xi, zi, wall){

    const iCoords = rotateLineSegment(
        xi,
        zi + agentLength + RADIUS,
        xi,
        zi - agentLength - RADIUS,
        sceneEntities[i].agent.rotation.z
    );



    // Agent A
    const a = {
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      radius: RADIUS,
    };
    // Wall B
    const b = {
      tip: new THREE.Vector3(wall.tip.x, 0, wall.tip.z),
      base: new THREE.Vector3(wall.base.x, 0, wall.base.z),
      radius: RADIUS,
    };


    // capsule A:
    const a_Normal = a.tip.clone().sub(a.base.clone()).normalize();
    const a_LineEndOffset = a_Normal.clone().multiplyScalar(a.radius);
    const a_A = a.base.clone().add(a_LineEndOffset);
    const a_B = a.tip.clone().sub(a_LineEndOffset);

    // capsule B:
    const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
    const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
    const b_A = b.base.clone().add(b_LineEndOffset);
    const b_B = b.tip.clone().sub(b_LineEndOffset);

    // vectors between line endpoints:
    const v0 = b_A.clone().sub(a_A);
    const v1 = b_B.clone().sub(a_A);
    const v2 = b_A.clone().sub(a_B);
    const v3 = b_B.clone().sub(a_B);

    // squared distances:
    const d0 = v0.clone().dot(v0);
    const d1 = v1.clone().dot(v1);
    const d2 = v2.clone().dot(v2);
    const d3 = v3.clone().dot(v3);

    // select best potential endpoint on capsule A:
    let bestA;
    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
      bestA = a_B;
    } else {
      bestA = a_A;
    }

    // select point on capsule B line segment nearest to best potential endpoint on A capsule:
    const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);

    // now do the same for capsule A segment:
    bestA = ClosestPointOnLineSegment(a_A, a_B, bestB);

    return [bestA, bestB, a, b]
  }



  function getBestPoint(xi, zi, xj, zj){

    const iCoords = rotateLineSegment(
        xi,
        zi + agentLength + RADIUS,
        xi,
        zi - agentLength - RADIUS,
        sceneEntities[i].agent.rotation.z
    );

    const jCoords = rotateLineSegment(
        xj,
        zj + agentLength + RADIUS,
        xj,
        zj - agentLength - RADIUS,
        sceneEntities[j].agent.rotation.z
    );


    // Agent A
    const a = {
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };
    // Agent B
    const b = {
      tip: new THREE.Vector3(jCoords[0], 0, jCoords[1]),
      base: new THREE.Vector3(jCoords[2], 0, jCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };


    // capsule A:
    const a_Normal = a.tip.clone().sub(a.base.clone()).normalize();
    const a_LineEndOffset = a_Normal.clone().multiplyScalar(a.radius);
    const a_A = a.base.clone().add(a_LineEndOffset);
    const a_B = a.tip.clone().sub(a_LineEndOffset);
    a.real_tip = a_B;
    a.real_base = a_A;

    // capsule B:
    const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
    const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
    const b_A = b.base.clone().add(b_LineEndOffset);
    const b_B = b.tip.clone().sub(b_LineEndOffset);
    b.real_tip = b_B;
    b.real_base = b_A;

    // vectors between line endpoints:
    const v0 = b_A.clone().sub(a_A);
    const v1 = b_B.clone().sub(a_A);
    const v2 = b_A.clone().sub(a_B);
    const v3 = b_B.clone().sub(a_B);

    // squared distances:
    const d0 = v0.clone().dot(v0);
    const d1 = v1.clone().dot(v1);
    const d2 = v2.clone().dot(v2);
    const d3 = v3.clone().dot(v3);

    // select best potential endpoint on capsule A:
    let bestA;
    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
      bestA = a_B;
    } else {
      bestA = a_A;
    }

    // select point on capsule B line segment nearest to best potential endpoint on A capsule:
    const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);

    // now do the same for capsule A segment:
    bestA = ClosestPointOnLineSegment(a_A, a_B, bestB);

    return [bestA, bestB, a, b]
  }

  function dotProduct(vector1, vector2) {
    let result = 0;
    for (let i = 0; i < vector1.length; i++) {
      result += vector1[i] * vector2[i];
    }
    return result;
  }

  function areCollinear(vector1, vector2) {
    // Ensure vectors are of the same dimension
    if (vector1.length !== vector2.length) {
      return false;
    }

    // Find the ratio of the first non-zero pair of elements
    let ratio;
    for (let i = 0; i < vector1.length; i++) {
      if (vector1[i] !== 0 && vector2[i] !== 0) {
        ratio = vector1[i] / vector2[i];
        break;
      }
    }

    // Check if all corresponding elements are in the same ratio
    for (let i = 0; i < vector1.length; i++) {
      // Handle division by zero cases
      if (vector1[i] === 0 && vector2[i] !== 0 || vector1[i] !== 0 && vector2[i] === 0) {
        return false;
      }

      // Check the ratio
      if (vector1[i] !== 0 && vector2[i] !== 0) {
        if (vector1[i] / vector2[i] !== ratio) {
          return false;
        }
      }
    }

    return true;
  }

  function signNoP(n){
    if (n >= 0){
      return 1;
    }else {
      return -1;
    }

  }


  /*  -----------------------  */

  agentVelocityPlanner();

  sceneEntities.forEach(function (item) {
    item.px = item.x + timestep * item.vx;
    item.pz = item.z + timestep * item.vz;
    item.py = item.y + timestep * item.vy;
  });






  let pbdIters = 0;
  let isColliding;
  var agent_a,
    agent_b,
    desDistance,
    i,
    j,
    idx = 0;

  while (pbdIters < ITERNUM) {


    // clean previous accumulated gradient
    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {

        sceneEntities[i].grad.x = 0;
        sceneEntities[i].grad.z = 0;
        sceneEntities[j].grad.x = 0;
        sceneEntities[j].grad.z = 0;

        sceneEntities[i].grad.dx = 0;
        sceneEntities[i].grad.dz = 0;
        sceneEntities[j].grad.dx = 0;
        sceneEntities[j].grad.dz = 0;

        j += 1;
      }
      i += 1;
    }



/*
    // wall collision (based on short range)
    i=0;
    while(i<sceneEntities.length)
    {
      j=0;
      while(j<customParams.wallData.length)
      {
        let [p_bestA, w_bestB, p_agent_i,p_agent_j] = getBestPointWithWall(sceneEntities[i].px, sceneEntities[i].pz, customParams.wallData[j]);

        let dist_wall_to_agent = distance(p_bestA.x, p_bestA.z, w_bestB.x, w_bestB.z);
        // let dist_cur_to_path_middle = distance(sceneEntities[i].goal_x, sceneEntities[i].goal_z, sceneEntities[i].px, sceneEntities[i].pz);

        // let dist_cur_to_path_middle = distance(-4, 4, sceneEntities[i].px, sceneEntities[i].pz);
        let dist_cur_to_path_middle = distance(sceneEntities[i].goal_x, sceneEntities[i].goal_z, sceneEntities[i].px, sceneEntities[i].pz);

        // let angle = calculateAngle(-4, 4, -4, -140, sceneEntities[i].vx, sceneEntities[i].vz);
        // let angle = calculateAngle(sceneEntities[i].vx, sceneEntities[i].vz, -4, -140, -4, 4);
        let angle = calculateAngle( -4, -140, -4, 4, sceneEntities[i].vx, sceneEntities[i].vz);

        // Calculate position vector connecting the two points
        let positionVector = { x: -4 - (-4), y: 4 - (-140) }
        // let positionVector = { x: -4 - (sceneEntities[i].goal_x), y: 4 - (-40) }
        let velocityVector = { x: sceneEntities[i].vx, y: sceneEntities[i].vz }; // Velocity vector of the agent

        // Calculate the angle between the position vector and the velocity vector
        let angle_2 = angleBetweenVectors_2(positionVector, velocityVector);

        let x_dist_diff = Math.abs(-4 - sceneEntities[i].px);
        let x_dist_diff_with_goal = Math.abs(sceneEntities[i].px - sceneEntities[i].goal_x);

        // if(dist_cur_to_path_middle < 13 && 180-angle_2 > 0.1)
        if(dist_cur_to_path_middle < 13 )
        {
          const angleRad = 0.5 * (Math.PI/180);
          // console.log("angleRad: ", angleRad);

          // sceneEntities[i].goal_x =  sceneEntities[i].goal_x - 2 * sceneEntities[i].goal_x * Math.cos(angleRad);
          // sceneEntities[i].goal_z =  sceneEntities[i].goal_z - 2 * sceneEntities[i].goal_z * Math.sin(angleRad);

          // sceneEntities[i].goal_x =  sceneEntities[i].goal_x - 4 * Math.cos(angleRad);
          sceneEntities[i].goal_z =  sceneEntities[i].goal_z - 4 * Math.sin(angleRad);
        }

        let penetration_normal = p_bestA.clone().sub(w_bestB);
        const len = penetration_normal.length();
        penetration_normal.divideScalar(len); // normalize
        const penetration_depth = sceneEntities[i].radius + Math.sqrt(2) * sceneEntities[i].radius - len;
        const intersects = penetration_depth > 0;
        if (intersects) {
          sceneEntities[i].colliding = true;

          sceneEntities[i].px += penetration_normal.x * 1 * penetration_depth;
          sceneEntities[i].pz += penetration_normal.z * 1 * penetration_depth;

        }

        j+=1;
      }
      i+=1
    }
*/



//=========================================== our Long-range call started ======================================================================

    //Capsule to Capsule long-range collision avoidance.
    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {
        let [bestA, bestB, agent_i, agent_j] = getBestPoint(sceneEntities[i].x, sceneEntities[i].z, sceneEntities[j].x, sceneEntities[j].z);
        let [p_bestA, p_bestB, p_agent_i,p_agent_j] = getBestPoint(sceneEntities[i].px, sceneEntities[i].pz, sceneEntities[j].px, sceneEntities[j].pz);
        // // ttc in long range collision paper

        let [delta_correction_i, delta_correction_j, grad_i, grad_j, s] = longRangeConstraintCapsule(
            bestA, bestB,
            p_bestA, p_bestB,
            sceneEntities[i].agent.rotation.z, sceneEntities[j].agent.rotation.z,
            agent_i, agent_j,
            sceneEntities[i], sceneEntities[j],
            i, j
        );

        sceneEntities[i].px += delta_correction_i.x;
        sceneEntities[i].pz += delta_correction_i.y;
        sceneEntities[j].px += delta_correction_j.x;
        sceneEntities[j].pz += delta_correction_j.y;

// -------------- Start --------------------- for calculating normal to capsule body --- Capsule A ----------------------
    // For capsule A
    let iCoords = rotateLineSegment(
      sceneEntities[i].x,
      sceneEntities[i].z + agentLength + RADIUS,
      sceneEntities[i].x,
      sceneEntities[i].z - agentLength - RADIUS,
      sceneEntities[i].agent.rotation.z
    );

    // Agent A
    const aa = {
    tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
    base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
    radius: RADIUS,
    real_tip: null,
    real_base: null
    };

    // capsule A:
    let a_Normal = aa.tip.clone().sub(aa.base.clone()).normalize();
    const a_LineEndOffset = a_Normal.clone().multiplyScalar(aa.radius);
    const a_A = aa.base.clone().add(a_LineEndOffset);
    const a_B = aa.tip.clone().sub(a_LineEndOffset);
    aa.real_tip = a_B;
    aa.real_base = a_A;

    // Calculate the slope of the line
    const dx = aa.base.x - aa.tip.x;
    const dz = aa.base.z - aa.tip.z  ;

    // Get the length of the line segment
    let leng = Math.sqrt(dx*dx + dz*dz);

    // Generate a normal vector (-dy, dx) 
    // This will be perpendicular to the line
    let nx = -dz / leng;  
    let nz = dx / leng;

    //------------end-------------------------------------------------

    // Normal vector to casule body
    let normal_to_capsule = new THREE.Vector3(nx, 0, nz);
    sceneEntities[i].normal_to_capsule = normal_to_capsule;
  
    var currentPosition = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
    var predicted_position = new THREE.Vector3(sceneEntities[i].px, 0, sceneEntities[i].pz);
    var displacementVector = predicted_position.clone().sub(currentPosition);
    var directionOfMovement = displacementVector.clone().normalize();
    var directionVector = new THREE.Vector3(directionOfMovement.x, 0, directionOfMovement.z);

    let angle_6 = angleBetweenVectors_2(directionVector, normal_to_capsule);

// =============== restrict big orientation change for Capsule A ==== Start =======================  
    if(angle_6 > angle_threshold && sceneEntities[i].x < 20 )
    {
      // console.log("before, angle_6: ", angle_6);
      let normal_to_capsule_vec = new THREE.Vector3(normal_to_capsule.x, 0, normal_to_capsule.z);
      let normal_to_capsule_vec_normalized = normal_to_capsule_vec.clone().normalize();

      let theta = angle_threshold * (Math.PI/180);
      let w_r = getNewVectorDeviatedFromAnother(normal_to_capsule_vec_normalized, theta);
      let w_r_direc = new THREE.Vector3(w_r.x, 0, w_r.z);
      w_r_direc = w_r_direc.clone().normalize();

      let new_point_right =  findNewPositionByProjecttion(predicted_position, currentPosition, w_r_direc);
     
      //checking whether projected point is correct or not by checking angle with new position vector and capsule normal.
      // var result2 = new THREE.Vector3();
      // result2.subVectors(currentPosition, new_point_right);
      // result2.clone().normalize();
      // let angl2 = angleBetweenVectors_2(result2, normal_to_capsule_vec_normalized);
      // console.log("result2: ", result2);
      // console.log("angl2: ", new_point_right, "");

      let theta_left = -angle_threshold * (Math.PI/180);
      let w_l = getNewVectorDeviatedFromAnother(normal_to_capsule_vec_normalized, theta_left);
      let w_l_direc = new THREE.Vector3(w_l.x, 0, w_l.z);
      w_l_direc = w_l_direc.clone().normalize();

      let new_point_left =  findNewPositionByProjecttion(predicted_position, currentPosition, w_l_direc);

      // var result = new THREE.Vector3();
      // result.subVectors(currentPosition, new_point_left);
      // let angl = angleBetweenVectors_2(result, normal_to_capsule);

      let temp_dist_right = distance(predicted_position.x, predicted_position.z, new_point_right.x, new_point_right.z);
      let temp_dist_left = distance(predicted_position.x, predicted_position.z, new_point_left.x, new_point_left.z);

      if(temp_dist_left >= temp_dist_right)
      {
        sceneEntities[i].px = new_point_right.x;
        sceneEntities[i].pz = new_point_right.z;
      }
      else{
          sceneEntities[i].px = new_point_left.x;
          sceneEntities[i].pz = new_point_left.z;
        }


      //checking current angle (after fix)
    // var currentPosition = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
    // var predicted_position = new THREE.Vector3(sceneEntities[i].px, 0, sceneEntities[i].pz);
    // var displacementVector = predicted_position.clone().sub(currentPosition);
    // var directionOfMovement = displacementVector.clone().normalize();
    // var directionVector = new THREE.Vector3(directionOfMovement.x, 0, directionOfMovement.z);

    // let angle_after = angleBetweenVectors_2(directionVector, normal_to_capsule);
    // console.log("angle_6 after : ", angle_after, "\n "); 
    }

// =============== restrict big orientation change for Capsule A ==== End =======================    

    
    // ------------------------------------------------------------------------------------------

    // For capsule B
    let jCoords = rotateLineSegment(
      sceneEntities[j].x,
      sceneEntities[j].z + agentLength + RADIUS,
      sceneEntities[j].x,
      sceneEntities[j].z - agentLength - RADIUS,
      sceneEntities[j].agent.rotation.z
    );

    // Agent B
    const bb = {
      tip: new THREE.Vector3(jCoords[0], 0, jCoords[1]),
      base: new THREE.Vector3(jCoords[2], 0, jCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
      };

    // capsule B
    let b_Normal = bb.tip.clone().sub(bb.base.clone()).normalize();
    const b_LineEndOffset = b_Normal.clone().multiplyScalar(bb.radius);
    const b_B = bb.base.clone().add(b_LineEndOffset);
    const a_B_2 = bb.tip.clone().sub(b_LineEndOffset);
    bb.real_tip = a_B_2;
    bb.real_base = b_B;

    // Calculate the slope of the line
    const dx_b = bb.base.x - bb.tip.x;
    const dz_b = bb.base.z - bb.tip.z  ;

    // Get the length of the line segment
    let leng_b = Math.sqrt(dx_b*dx_b + dz_b*dz_b);

    // Generate a normal vector (-dy, dx) 
    // This will be perpendicular to the line
    let nx_b = -dz_b / leng_b;  
    let nz_b = dx_b / leng_b;

    // Normal vector to capsule body
    let normal_to_capsule_b = new THREE.Vector3(nx_b, 0, nz_b);
    normal_to_capsule_b = normal_to_capsule_b.clone().normalize();

    sceneEntities[j].normal_to_capsule = normal_to_capsule_b;

    var currentPosition = new THREE.Vector3(sceneEntities[j].x, 0, sceneEntities[j].z);
    var predicted_position = new THREE.Vector3(sceneEntities[j].px, 0, sceneEntities[j].pz);
    var displacementVector = predicted_position.clone().sub(currentPosition);
    var directionOfMovement = displacementVector.clone().normalize();
    var directionVector = new THREE.Vector3(directionOfMovement.x, 0, directionOfMovement.z);

    let angle_4_b = angleBetweenVectors_2(directionVector, normal_to_capsule_b);
    // console.log("angle_4_b: ", angle_4_b);


// =============== restrict big orientation change for Capsule B ==== Start =======================

  if(angle_4_b > angle_threshold)
  {
    // console.log("before, angle_4_b: ", angle_4_b);
    let normal_to_capsule_vec_b = new THREE.Vector3(normal_to_capsule_b.x, 0, normal_to_capsule_b.z);
    let normal_to_capsule_vec_normalized_b = normal_to_capsule_vec_b.clone().normalize();

    let theta_right = angle_threshold * (Math.PI/180);
    let w_r_b = getNewVectorDeviatedFromAnother(normal_to_capsule_vec_normalized_b, theta_right);
    let w_r_direc_b = new THREE.Vector3(w_r_b.x, 0, w_r_b.z);
    w_r_direc_b = w_r_direc_b.clone().normalize();

    let new_point_right =  findNewPositionByProjecttion(predicted_position, currentPosition, w_r_direc_b);

    //  check whether angle after projecting predicted position over direction vector right or not.
    //  var result = new THREE.Vector3();
    //  result.subVectors(new_point_right, currentPosition);
    //  let angl = angleBetweenVectors_2(normal_to_capsule_vec_normalized_b, result);
    //  console.log("angl: ", angl);

    let theta_left = -angle_threshold * (Math.PI/180);
    let w_l = getNewVectorDeviatedFromAnother(normal_to_capsule_vec_normalized_b, theta_left);
    let w_l_direc = new THREE.Vector3(w_l.x, 0, w_l.z);
    w_l_direc = w_l_direc.clone().normalize();

    let new_point_left =  findNewPositionByProjecttion(predicted_position, currentPosition, w_l_direc);

    let temp_dist_right = distance(predicted_position.x, predicted_position.z, new_point_right.x, new_point_right.z);
    let temp_dist_left = distance(predicted_position.x, predicted_position.z, new_point_left.x, new_point_left.z);

    if(temp_dist_left > temp_dist_right)
    {
      sceneEntities[j].px = new_point_right.x;
      sceneEntities[j].pz = new_point_right.z;
    }
    else{
        sceneEntities[j].px = new_point_left.x;
        sceneEntities[j].pz = new_point_left.z;
      }


  //checking current angle (after fix)
// var currentPosition = new THREE.Vector3(sceneEntities[j].x, 0, sceneEntities[j].z);
// var predicted_position = new THREE.Vector3(sceneEntities[j].px, 0, sceneEntities[j].pz);
// var displacementVector = predicted_position.clone().sub(currentPosition);
// var directionOfMovement = displacementVector.clone().normalize();
// var directionVector = new THREE.Vector3(directionOfMovement.x, 0, directionOfMovement.z);

// let angle_after = angleBetweenVectors_2(directionVector, normal_to_capsule_b);
// console.log("angle_4_b after : ", angle_after, "\n ");
}

// =============== restrict big orientation change for Capsule B ===== End =======================

// ------------------ END ----------------- for calculating normal to capsule body ----------------------


        // for utilities
        sceneEntities[i].grad.x += grad_i[0];
        sceneEntities[i].grad.z += grad_i[1];
        sceneEntities[j].grad.x += grad_j[0];
        sceneEntities[j].grad.z += grad_j[1];

        sceneEntities[i].grad.s = s;
        sceneEntities[j].grad.s = s;

        sceneEntities[i].grad.dx += delta_correction_i.x;
        sceneEntities[i].grad.dz += delta_correction_i.y;
        sceneEntities[j].grad.dx += delta_correction_j.x;
        sceneEntities[j].grad.dz += delta_correction_j.y;

        customParams.best[i][j] = [bestA, bestB]
        customParams.best[j][i] = [bestB, bestA]

        // short range collision
        // didn't correct position in real time
        // let penetration_normal = bestA.clone().sub(bestB);
        // const len = penetration_normal.length();
        // penetration_normal.divideScalar(len); // normalize
        // const penetration_depth = sceneEntities[i].radius + sceneEntities[j].radius - len;
        // const intersects = penetration_depth > 0;
        // if (intersects) {
        //   sceneEntities[i].colliding = true;
        //   sceneEntities[j].colliding = true;
        //
        //   sceneEntities[i].px += penetration_normal.x * 0.5 * penetration_depth;
        //   sceneEntities[i].pz += penetration_normal.y * 0.5 * penetration_depth;
        //
        //   sceneEntities[j].px +=
        //       -1 * penetration_normal.x * 0.5 * penetration_depth;
        //   sceneEntities[j].pz +=
        //       -1 * penetration_normal.y * 0.5 * penetration_depth;
        // }

        j += 1;
      }

      i += 1;
    }
//======================================== our Long-range call ended ===============================================================




  
    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {

        let [bestA, bestB, agent_i, agent_j] = getBestPoint(sceneEntities[i].x, sceneEntities[i].z, sceneEntities[j].x, sceneEntities[j].z);
        let [p_bestA, p_bestB, p_agent_i,p_agent_j] = getBestPoint(sceneEntities[i].px, sceneEntities[i].pz, sceneEntities[j].px, sceneEntities[j].pz);
        
        let [delta_correction_i, delta_correction_j, grad_i, grad_j, s] = collisionConstraint_Capsule(bestA, bestB, p_bestA, p_bestB);
        j += 1;
      }
      i += 1;
    }


    pbdIters += 1;
  }


  sceneEntities.forEach(function (item) {

    // const dx = item.goal_x - item.x;
    // const dz = item.goal_z - item.z;

    const dx = item.px - item.x;
    const dz = item.pz - item.z;

    item.agent.rotation.z = Math.atan2(dz, dx); 

    item.vx = (item.px - item.x) / timestep;
    item.vz = (item.pz - item.z) / timestep;
    item.vy = (item.py - item.y) / timestep;

    item.x = item.px;
    item.z = item.pz;
    item.y = item.py;
  });

}
