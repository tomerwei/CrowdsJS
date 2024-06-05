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


// Function to get a vector direction deviated by a certain angle from a known normalized vector
function deviateVectorByAngle(knownVector, angle) {

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


function findNewPositionByProjection(predict_point, cur_point, direction_vector, angle_btwn_cur_direc_and_capsule_body_normal) {

  let position_vector = {};

  position_vector = new THREE.Vector3().subVectors(predict_point, cur_point);

  // if(angle_btwn_cur_direc_and_capsule_body_normal < 90)
  // {
  // // get position vector from current and precidted position
  //  position_vector = new THREE.Vector3().subVectors(predict_point, cur_point);
  // }else{
  //   position_vector = new THREE.Vector3().subVectors(cur_point, predict_point);
  // }

  // Project position_vector onto w_r
  const proj = position_vector.clone().projectOnVector(direction_vector); 
  
  // Calculate new predicted positionon on w_r vector direction.
  const p_new = cur_point.clone().add(proj);

  return p_new;  
}


let closest_in_left = new THREE.Vector3();
let closest_in_right = new THREE.Vector3();

// let closest_in_left_agent = new THREE.Vector3();
// let closest_in_right_agent = new THREE.Vector3();
// let closest_in_left_wall = new THREE.Vector3();
// let closest_in_right_wall = new THREE.Vector3();

let best_in_left_agent = new THREE.Vector3();
let best_in_right_agent = new THREE.Vector3();
let best_in_left_wall = new THREE.Vector3();
let best_in_right_wall = new THREE.Vector3();

let new_point = new THREE.Vector3();
let normalized_velocity = {};

export function step(RADIUS, sceneEntities, world, scene, customParams = {}) {
  const AGENTSIZE = RADIUS * 2;
  const epsilon = 0.0001;
  const timestep = 0.02;
  const ITERNUM = 1; // 3
  const agentLength = RADIUS;

  let C_TAU_MAX = 20;
  let C_TAO0 = 250; 

  const C_LONG_RANGE_STIFF = 0.15;  
  const MAX_DELTA = 0.01;

  let angleThresholdBtwnDirectionAndNormalInDeg = 0.08; 

  let dist_tip_to_base  = 0;

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


function getCapsuleBodyNormal(agent, agentLength, RADIUS, current_rotation) {

    let iCoords = rotateLineSegment(
        agent.x,
        agent.z + agentLength + RADIUS,
        agent.x,
        agent.z - agentLength - RADIUS,
        current_rotation
      );

    if(customParams.scenario==='bottleneck')
    {
       iCoords = rotateLineSegment(
        agent.z,
        agent.x + agentLength + RADIUS,
        agent.z,
        agent.x - agentLength - RADIUS,
        current_rotation
      );
    }
  
    const aa = {
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      };
      
      // console.log("dist: ", distance(aa.tip, aa.base));

      // Calculate the slope of the line
      const dx = aa.base.x - aa.tip.x;
      const dz = aa.base.z - aa.tip.z  ;
  
      let leng = Math.sqrt(dx*dx + dz*dz);
        
      // Generate a normal vector (-dy, dx), perpendicular to the line
      let nx = -dz / leng;  
      let nz = dx / leng;
  
      let normal_to_capsule = new THREE.Vector3(nx, 0, nz);
      agent.normal_to_capsule = normal_to_capsule;

      
      let angle_capsule_normal_and_vel = angleBetweenVectors_2(agent.normal_to_capsule, agent.normal_to_capsule_prev);

      if(angle_capsule_normal_and_vel > 20){    // forcing normal vector to capsule body to be in the capsule's facing direction. 
        agent.normal_to_capsule = agent.normal_to_capsule_prev;
      }
      agent.normal_to_capsule_prev = agent.normal_to_capsule;
      
  
    return agent.normal_to_capsule;
  }


  /*  -----------------------  */
  /*  TODO modify lines below  */
  /*  -----------------------  */


  function collisionConstraint_Capsule(best_i, best_j, p_best_i, p_best_j) {
    // let stif = 0.10;
    // let stif = 0.60;
    // let stif = 0.85;
    let stif = 1.0;

    const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

    const agentDist = agentCentroidDist - AGENTSIZE;
    const dir_x = (p_best_j.x - p_best_i.x) / agentCentroidDist;
    const dir_z = (p_best_j.z - p_best_i.z) / agentCentroidDist;
    const agent_i_scaler = (1 / (1 + 1)) * agentDist;
    const agent_j_scaler = (1 / (1 + 1)) * agentDist;

    if (agentDist < 0) {
      // console.log("capsules collided!! ", "");
      sceneEntities[i].px += stif *  agent_i_scaler * dir_x;
      sceneEntities[i].pz += stif * agent_i_scaler * dir_z;
      sceneEntities[j].px += stif * -agent_j_scaler * dir_x;
      sceneEntities[j].pz += stif * -agent_j_scaler * dir_z;
        
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

      // best points distance difference
      let projectedPoint_j = PointOnLineSegment(agent_j.real_base, agent_j.real_tip, best_i);
      let bestPointDiff = distance(projectedPoint_j.x, projectedPoint_j.z, best_j.x, best_j.z);

      // if facing direction on the same line AND the best points are exactly facing with each other
      // adding gradient value
  
      // if (bestPointDiff <= 4.0 &&  (facingDiff < 0.015 || facingDiff - Math.PI) <0.015 ){  //used this
      if (bestPointDiff <= 6.0 &&  (facingDiff < 0.015 || facingDiff - Math.PI) <0.015 ){
        
        grad_y_i = signNoP(grad_y_i) * (Math.abs(grad_x_i)/1.2);   //used this
        grad_y_j = -grad_y_i;
      }
    
      const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
      s = stiff * tao_sq / (0.5 * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + 0.5 * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed

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



  function rotationConstraint(capsule_entity)
  {
    let current_rotation = capsule_entity.agent.rotation.z;
  
    let NormalVectorToCapsuleBody = getCapsuleBodyNormal(capsule_entity, agentLength, RADIUS, current_rotation);
    let NormalVectorToCapsuleBodyNormalized = NormalVectorToCapsuleBody.clone().normalize();
  
    var currentPosition = new THREE.Vector3(capsule_entity.x, 0, capsule_entity.z);
    var predicted_position = new THREE.Vector3(capsule_entity.px, 0, capsule_entity.pz);
    var directionVector = predicted_position.clone().sub(currentPosition);    
    var directionVectorNormalized = directionVector.clone().normalize();
  
    let angleOfDirectionAndCapsuleBodyNormal = angleBetweenVectors_2(directionVectorNormalized, NormalVectorToCapsuleBody);
    
    //check if the angle btwn capsule body normal vector is greater than a threhold angle
    if(angleOfDirectionAndCapsuleBodyNormal > angleThresholdBtwnDirectionAndNormalInDeg )
    {  
      //calculating projected point on right side
      let angleOfDirectionAndCapsuleBodyNormal_Radians_right = angleThresholdBtwnDirectionAndNormalInDeg * (Math.PI/180);
      let rightSideDirectionVector_temp = deviateVectorByAngle(NormalVectorToCapsuleBodyNormalized, angleOfDirectionAndCapsuleBodyNormal_Radians_right);
      let rightSideDirectionVector = new THREE.Vector3(rightSideDirectionVector_temp.x, 0, rightSideDirectionVector_temp.z);
      rightSideDirectionVector = rightSideDirectionVector.clone().normalize();

      let projectedPointOnRightSide =  findNewPositionByProjection(predicted_position, currentPosition, rightSideDirectionVector, angleOfDirectionAndCapsuleBodyNormal);

      //calculating projected point on left side
      let angleOfDirectionAndCapsuleBodyNormal_Radians_Left= -angleThresholdBtwnDirectionAndNormalInDeg * (Math.PI/180);
      let leftSideDirectionVector_temp = deviateVectorByAngle(NormalVectorToCapsuleBodyNormalized, angleOfDirectionAndCapsuleBodyNormal_Radians_Left);
      let leftSideDirectionVector = new THREE.Vector3(leftSideDirectionVector_temp.x, 0, leftSideDirectionVector_temp.z);
      leftSideDirectionVector = leftSideDirectionVector.clone().normalize();

      let projectedPointOnLeftSide = findNewPositionByProjection(predicted_position, currentPosition, leftSideDirectionVector, angleOfDirectionAndCapsuleBodyNormal);

      // calculating distance btwn new predicted positions with old predicted position in left and right side of the capsule normal vector
      let distPredictedAndProjectedRightSide = distance(predicted_position.x, predicted_position.z, projectedPointOnRightSide.x, projectedPointOnRightSide.z);
      let distPredictedAndProjectedLeftSide = distance(predicted_position.x, predicted_position.z, projectedPointOnLeftSide.x, projectedPointOnLeftSide.z);
      
      //update agent position with whichever projected point has the shortest distance with old predicted position
      if(distPredictedAndProjectedLeftSide >= distPredictedAndProjectedRightSide )
      {
        capsule_entity.px = projectedPointOnRightSide.x;
        capsule_entity.pz = projectedPointOnRightSide.z;
      }
      else{
        capsule_entity.px = projectedPointOnLeftSide.x;
        capsule_entity.pz = projectedPointOnLeftSide.z;
      }
    }
  
    return { x: capsule_entity.px, y: 0, z: capsule_entity.pz };
  }


  function findLeftPointWithCertainDistance(capsule_entity, DeviationAngle)
  {
    // let angleTesting = 30;
    let distanceRadius = 2;
    let current_rotation = capsule_entity.agent.rotation.z;
  
    let NormalVectorToCapsuleBody = getCapsuleBodyNormal(capsule_entity, agentLength, RADIUS, current_rotation);
    let NormalVectorToCapsuleBodyNormalized = NormalVectorToCapsuleBody.clone().normalize();
  
    //calculating projected point on right side
    let angleOfDirectionAndCapsuleBodyNormal_Radians_right = DeviationAngle * (Math.PI/180);
    let rightSideDirectionVector_temp = deviateVectorByAngle(NormalVectorToCapsuleBodyNormalized, angleOfDirectionAndCapsuleBodyNormal_Radians_right);
    let rightSideDirectionVector = new THREE.Vector3(rightSideDirectionVector_temp.x, 0, rightSideDirectionVector_temp.z);
    rightSideDirectionVector = rightSideDirectionVector.clone().normalize();

    let newX = capsule_entity.x + rightSideDirectionVector.x * distanceRadius;
    let newY = capsule_entity.y + rightSideDirectionVector.y * distanceRadius;
    let newZ = capsule_entity.z + rightSideDirectionVector.z * distanceRadius;

    // let newPoint = new THREE.Vector3(newX, newY, newZ);

    return { x: newX, y: 0, z: newZ };
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

  
    // Calculate the slope of the line
    const dx = a.base.x - a.tip.x;
    const dz = a.base.z - a.tip.z  ;
    
    let leng = Math.sqrt(dx*dx + dz*dz);
        
    dist_tip_to_base = leng;

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



  function find_best_point_On_right_or_Left_of_Velocity(capsule_entity, capsule_or_wall_entity, agents_pair_type)
  {
      let bestA_agent = new THREE.Vector3();
      let bestB_agent = new THREE.Vector3();
      let bestA2 = new THREE.Vector3();
      let bestB_agent_or_wall = new THREE.Vector3();

      let cur_pos = new THREE.Vector3(capsule_entity.x, 0, capsule_entity.z);
      normalized_velocity = new THREE.Vector3(sceneEntities[i].vx, 0, sceneEntities[i].vz).normalize();
      let threshold_Dist = 0.5;
      
      // Calculate the new position
      const displacement = normalized_velocity.clone().multiplyScalar(threshold_Dist);
      new_point = cur_pos.clone().add(displacement);
      
      if(agents_pair_type == "agents_only")
      {
        [bestA_agent, bestB_agent_or_wall, , ] = getBestPoint(new_point.x, new_point.z, capsule_or_wall_entity.x, capsule_or_wall_entity.z);
       }else if(agents_pair_type == "agents_and_walls")
      {
        [bestA_agent, bestB_agent_or_wall, , ] = getBestPointWithWall(new_point.x, new_point.z, capsule_or_wall_entity);
       }
       
      let point_on_obs = new THREE.Vector3(bestB_agent_or_wall.x, 0, bestB_agent_or_wall.z);
      
      // Create vector from capsule position to the point
      const vectorToPoint = new THREE.Vector3().subVectors(point_on_obs, cur_pos);    
  
      // Rotate the velocity vector by 90 degrees in the 2D plane to get the perpendicular vector
      const perpendicularVector = new THREE.Vector3(-normalized_velocity.z, 0, normalized_velocity.x);
  
      // Compute the dot product
      const dotProduct = perpendicularVector.dot(vectorToPoint);
  
      // Determine if the point is to the left or right of the velocity vector.
      const direction = dotProduct > 0 ? 'right' : 'left';
  
      if(direction == 'right')
      {
          // closest_in_right = bestB_agent_or_wall;    // point is on the right side
          best_in_right_wall = bestB_agent_or_wall;
      }
  
      if(direction == 'left')             
      {
        // closest_in_left = bestB_agent_or_wall;       // point is on the left side
        best_in_left_wall = bestB_agent_or_wall;
      }
  
  return [best_in_right_wall, best_in_left_wall]
  }


  function findClosestInLeftAndRightOfVelocity(new_point, best_in_left_wall, best_in_left_agent, best_in_right_wall, best_in_right_agent)
  {
  
      let cur_to_left_wall_dist = distance(new_point.x, new_point.z, best_in_left_wall.x, best_in_left_wall.z);
      let cur_to_left_agent_dist = distance(new_point.x, new_point.z, best_in_left_agent.x, best_in_left_agent.z);
  
      let cur_to_right_wall_dist = distance(new_point.x, new_point.z, best_in_right_wall.x, best_in_right_wall.z);
      let cur_to_right_agent_dist = distance(new_point.x, new_point.z, best_in_right_agent.x, best_in_right_agent.z);
  
      if(cur_to_left_wall_dist >= cur_to_left_agent_dist)
      {
          closest_in_left = best_in_left_agent;
      }else{
          {
            closest_in_left = best_in_left_wall;
          }
      }
  
      if(cur_to_right_wall_dist >= cur_to_right_agent_dist)
      {
          closest_in_right = best_in_right_agent;
        }else{
          {
            closest_in_right = best_in_right_wall;
          }
      }
  
  return [closest_in_right, closest_in_left]
  }
  
  
  function compute_Shortest_Perpendicular_dist(closest_in_left_or_right, new_point, normalized_velocity_2){
  
        // Calculate the vector from the line origin to the point in question
        const originToPoint = new THREE.Vector3().subVectors(closest_in_left_or_right, new_point);
  
        // Project this vector onto the direction vector of the line
        const projection = normalized_velocity_2.clone().multiplyScalar(originToPoint.dot(normalized_velocity_2));
  
        // Calculate the perpendicular vector from the point to the line
        const perpendicularVector = new THREE.Vector3().subVectors(originToPoint, projection);
  
        // The length of this perpendicular vector is the distance from the point to the line
        let perpendicularDistance_right_or_left = perpendicularVector.length();
  
  return perpendicularDistance_right_or_left;
  }

  function computerClearanceForAgent(normalized_velocity_2, closest_in_right, closest_in_left)
  {
  
      let perpendicularDistance_right = 0;
      let perpendicularDistance_left = 0;
  
      let closest_in_right_2 = new THREE.Vector3(closest_in_right.x, 0, closest_in_right.z);
      let closest_in_left_2 = new THREE.Vector3(closest_in_left.x, 0, closest_in_left.z);
  
      let wall_depth = customParams.wallData[1].depth;
  
      if(closest_in_right_2.length() != 0)
      {
    perpendicularDistance_right = compute_Shortest_Perpendicular_dist(closest_in_right_2, new_point, normalized_velocity_2);
      }
  
      if(closest_in_left_2.length() != 0)
      {
    perpendicularDistance_left = compute_Shortest_Perpendicular_dist(closest_in_left_2, new_point, normalized_velocity_2);
      }


    let clearance = perpendicularDistance_right + perpendicularDistance_left - 2 * wall_depth;
  
  return clearance;
  }


  function orientationConstraint(capsule, clearance, dis_btwn_agents)
  {

    let cross_width = (dist_tip_to_base/2) + RADIUS;

    // if( (dis_agent_to_walls < 2.5) && clearance - sceneEntities[i].radius >= 3 * sceneEntities[i].radius )
    // if( (dis_btwn_agents < 5) && clearance - sceneEntities[i].radius >= cross_width )   //was before
    if( (dis_btwn_agents < (cross_width + 2)) && (clearance - capsule.radius) >= (2 * cross_width) )
    {
      customParams.orientation=  'side_step';
    // }else if( (dis_agent_to_walls < 2.5) && (0 < clearance - sceneEntities[i].radius) && ((clearance - sceneEntities[i].radius) < (3 )) ){
    // }else if( (dis_btwn_agents < 5) && ( clearance - sceneEntities[i].radius > 0) && ((clearance - sceneEntities[i].radius) < (cross_width )) ){   //was before
    }else if( (dis_btwn_agents < (cross_width+2)) && ( clearance - capsule.radius > 0) && ((clearance - capsule.radius) < (2 * cross_width )) ){
      // let cosValue = ( (clearance - sceneEntities[i].radius) / (3));
      // let cosValue = ( (clearance - sceneEntities[i].radius) / (cross_width ));
      let cosValue = ( (clearance - capsule.radius) / (2 * cross_width ));
      const angleInRadians = Math.acos(cosValue);
      // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + angleInRadians;
      // sceneEntities[i].agent.rotation.z = angleInRadians;
      if( Math.abs(capsule.agent.rotation.z - angleInRadians) >= 0.08)
      {
        // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + angleInRadians/100; //was before
        capsule.agent.rotation.z = capsule.agent.rotation.z + angleInRadians/150;
      }
      
      customParams.orientation=  'side_step';

    // }else if( (dis_agent_to_walls < 2.5) && clearance - sceneEntities[i].radius <= 0){
    }else if( (dis_btwn_agents < (cross_width+2) ) && clearance - capsule.radius <= 0){
      // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + 1.5708;

    }else{
        customParams.orientation=  'front';
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


// agent to wall short-range collision constraint    
i=0;
   while(i<sceneEntities.length)
   {
     j=0;
     while(j<customParams.wallData.length)
     {
       let [p_bestA, w_bestB, p_agent_i,p_agent_j] = getBestPointWithWall(sceneEntities[i].px, sceneEntities[i].pz, customParams.wallData[j]);

       let penetration_normal = p_bestA.clone().sub(w_bestB);
       const len = penetration_normal.length();
       penetration_normal.divideScalar(len); // normalize
       const penetration_depth = sceneEntities[i].radius + 0.50 - len ;  //0.5 is the depth of the wall

       const intersects = penetration_depth > 0;
       if (intersects) {
         sceneEntities[i].colliding = true;
          sceneEntities[i].px += penetration_normal.x * 1.0 * penetration_depth;  
         sceneEntities[i].pz += penetration_normal.z * 1.0 * penetration_depth;
       }
       j+=1;
     }
     i+=1
   }





//Orientation constraint
//============== start ==  orientation constraint =========================== For one agent passing through walls scenario ========================================================================================================
i = 0;
let dis_btwn_agents = 0;
while (i < sceneEntities.length) {
  j = i ;

  let wall_depth = 0;
  let new_point = new THREE.Vector3();
  wall_depth = customParams.wallData[1].depth;

  let cur_pos = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
  normalized_velocity = new THREE.Vector3(sceneEntities[i].vx, 0, sceneEntities[i].vz).normalize();
  let threshold_Dist = 0.5;
  
  // Calculate the new position
  const displacement = normalized_velocity.clone().multiplyScalar(threshold_Dist);
  new_point = cur_pos.clone().add(displacement);
  
  while(j<customParams.wallData.length){
    
 /*   
    //----------------- start ------- Clearance computation ------ -------------------------------------------------
    let cur_pos = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
    normalized_velocity = new THREE.Vector3(sceneEntities[i].vx, 0, sceneEntities[i].vz).normalize();
    let threshold_Dist = 0.5;
    
    // Calculate the new position
    const displacement = normalized_velocity.clone().multiplyScalar(threshold_Dist);
    new_point = cur_pos.clone().add(displacement);
    
    let [bestA2, w_bestB2, agent_i, agent_j] = getBestPointWithWall(new_point.x, new_point.z, customParams.wallData[j]);  
    let point_on_obs = new THREE.Vector3(w_bestB2.x, 0, w_bestB2.z);
    
    // Create vector from capsule position to the point
    const vectorToPoint = new THREE.Vector3().subVectors(point_on_obs, cur_pos);    

    // Rotate the velocity vector by 90 degrees in the 2D plane to get the perpendicular vector
    const perpendicularVector = new THREE.Vector3(-normalized_velocity.z, 0, normalized_velocity.x);

    // Compute the dot product
    const dotProduct = perpendicularVector.dot(vectorToPoint);

    // Determine if the point is to the left or right of the velocity vector.
    const direction = dotProduct > 0 ? 'right' : 'left';

    if(direction == 'right')
    {
        // closest_in_right = w_bestB2;    // point is on the right side
        closest_in_right_wall = w_bestB2;
    }

    if(direction == 'left')             
    {
      // closest_in_left = w_bestB2;       // point is on the left side
      closest_in_left_wall = w_bestB2;
    }
*/
    [best_in_right_wall, best_in_left_wall] = find_best_point_On_right_or_Left_of_Velocity(sceneEntities[i], customParams.wallData[j], "agents_and_walls");

    j += 1;
  }



//----------------------------- START ------------------------------- to iterate over all agents ------------------------------------------------------------

if(customParams.scenario === 'narrow_hallwayOneAgent_Scenario'){
  closest_in_right = best_in_right_wall;
    closest_in_left = best_in_left_wall;
    console.log("testing");
}
if(customParams.scenario === 'narrow_hallwayTwoAgent_FaceToFace'){
    j = i+1 ;
    while (j < sceneEntities.length) {

      dis_btwn_agents = distance(sceneEntities[i].x, sceneEntities[i].z, sceneEntities[j].x, sceneEntities[j].z);
/*
      let cur_pos = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z);
      normalized_velocity = new THREE.Vector3(sceneEntities[i].vx, 0, sceneEntities[i].vz).normalize();
    
      let threshold_Dist = 0.5; 

      // Calculate the new position
      const displacement = normalized_velocity.clone().multiplyScalar(threshold_Dist);
      new_point = cur_pos.clone().add(displacement);

      let [bestA_agent, bestB_agent, agent_i2, agent_j2] = getBestPoint(new_point.x, new_point.z, sceneEntities[j].x, sceneEntities[j].z);   

      let other_point = new THREE.Vector3(bestB_agent.x, 0, bestB_agent.z);
      // Create vector from capsule position to the point
      const vectorToPoint = new THREE.Vector3().subVectors(other_point, cur_pos);        

      // Rotate the velocity vector by 90 degrees in the 2D plane to get the perpendicular vector
      const perpendicularVector = new THREE.Vector3(-normalized_velocity.z, 0, normalized_velocity.x);

      // Compute the dot product
      const dotProduct = perpendicularVector.dot(vectorToPoint);

      // Determine if the point is to the left or right
      const direction = dotProduct > 0 ? 'right' : 'left';
      
      if(direction == "right")
      {
          closest_in_right_agent = bestB_agent;
      }

      if(direction == "left"){
        closest_in_left_agent = bestB_agent;
    }
*/

    [best_in_right_agent, best_in_left_agent] = find_best_point_On_right_or_Left_of_Velocity(sceneEntities[i], sceneEntities[j], "agents_only");

      j += 1;
    }
//----------------------------- END ------------------------------- to iterate over all agents ------------------------------------------------------------
/*
    let cur_to_left_wall_dist = distance(new_point.x, new_point.z, closest_in_left_wall.x, closest_in_left_wall.z);
    let cur_to_left_agent_dist = distance(new_point.x, new_point.z, closest_in_left_agent.x, closest_in_left_agent.z);

    let cur_to_right_wall_dist = distance(new_point.x, new_point.z, closest_in_right_wall.x, closest_in_right_wall.z);
    let cur_to_right_agent_dist = distance(new_point.x, new_point.z, closest_in_right_agent.x, closest_in_right_agent.z);

    if(cur_to_left_wall_dist >= cur_to_left_agent_dist)
    {
        closest_in_left = closest_in_left_agent;
    }else{
        {
          closest_in_left = closest_in_left_wall;
        }
    }

    if(cur_to_right_wall_dist >= cur_to_right_agent_dist)
    {
        closest_in_right = closest_in_right_agent;
      }else{
        
        closest_in_right = closest_in_right_wall;
        
    }

    */

    [closest_in_right, closest_in_left] = findClosestInLeftAndRightOfVelocity(new_point, best_in_left_wall, best_in_left_agent, best_in_right_wall, best_in_right_agent)

    }

    let normalized_velocity_2 = new THREE.Vector3(normalized_velocity.x, 0, normalized_velocity.z);
    let clearance = computerClearanceForAgent(normalized_velocity_2, closest_in_right, closest_in_left);


/*    
    let perpendicularDistance_right = 0;
    let perpendicularDistance_left = 0;

    if(closest_in_right.length() != 0)
    {
      let closest_in_right_2 = new THREE.Vector3(closest_in_right.x, 0, closest_in_right.z);
      // let normalized_velocity_2 = new THREE.Vector3(normalized_velocity.x, 0, normalized_velocity.z);

      // Define the line using an origin point and a direction vector
      // const pointOnVelVector = new THREE.Vector3(sceneEntities[i].x, 0, sceneEntities[i].z); // Origin point on the line

      // Calculate the vector from the line origin to the point in question
      const originToPoint = new THREE.Vector3().subVectors(closest_in_right_2, new_point);

      // Project this vector onto the direction vector of the line
      const projection = normalized_velocity_2.clone().multiplyScalar(originToPoint.dot(normalized_velocity_2));

      // Calculate the perpendicular vector from the point to the line
      const perpendicularVector = new THREE.Vector3().subVectors(originToPoint, projection);

      // The length of this perpendicular vector is the distance from the point to the line
      perpendicularDistance_right = perpendicularVector.length();
    }

    if(closest_in_left.length() != 0)
    {
      let closest_in_left_2 = new THREE.Vector3(closest_in_left.x, 0, closest_in_left.z);
      // let normalized_velocity_2 = new THREE.Vector3(normalized_velocity.x, 0, normalized_velocity.z);

      // Calculate the vector from the line origin to the point in question
      const originToPoint = new THREE.Vector3().subVectors(closest_in_left_2, new_point);

      // Project this vector onto the direction vector of the line
      const projection = normalized_velocity_2.clone().multiplyScalar(originToPoint.dot(normalized_velocity_2));

      // Calculate the perpendicular vector from the point to the line
      const perpendicularVector = new THREE.Vector3().subVectors(originToPoint, projection);

      // The length of this perpendicular vector is the distance from the point to the line
      perpendicularDistance_left = perpendicularVector.length();
    }

    // let clearance = perpendicularDistance_right + perpendicularDistance_left - 1.4 * wall_depth;   // minus depth of the right and left walls to find the actual available clearance.
    let clearance = perpendicularDistance_right + perpendicularDistance_left - 2 * wall_depth;
    // let clearance = perpendicularDistance_right + perpendicularDistance_left; 

*/


/*
    let cross_width = (dist_tip_to_base/2) + RADIUS;

    // if( (dis_agent_to_walls < 2.5) && clearance - sceneEntities[i].radius >= 3 * sceneEntities[i].radius )
    // if( (dis_btwn_agents < 5) && clearance - sceneEntities[i].radius >= cross_width )   //was before
    if( (dis_btwn_agents < (cross_width+2)) && (clearance - sceneEntities[i].radius) >= (2 * cross_width) )
    {
      customParams.orientation=  'side_step';
    // }else if( (dis_agent_to_walls < 2.5) && (0 < clearance - sceneEntities[i].radius) && ((clearance - sceneEntities[i].radius) < (3 )) ){
    // }else if( (dis_btwn_agents < 5) && ( clearance - sceneEntities[i].radius > 0) && ((clearance - sceneEntities[i].radius) < (cross_width )) ){   //was before
    }else if( (dis_btwn_agents < (cross_width+2)) && ( clearance - sceneEntities[i].radius > 0) && ((clearance - sceneEntities[i].radius) < (2 * cross_width )) ){
      // let cosValue = ( (clearance - sceneEntities[i].radius) / (3));
      // let cosValue = ( (clearance - sceneEntities[i].radius) / (cross_width ));
      let cosValue = ( (clearance - sceneEntities[i].radius) / (2 * cross_width ));
      const angleInRadians = Math.acos(cosValue);
      // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + angleInRadians;
      // sceneEntities[i].agent.rotation.z = angleInRadians;
      if( Math.abs(sceneEntities[i].agent.rotation.z - angleInRadians) >= 0.08)
      {
        // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + angleInRadians/100; //was before
        sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + angleInRadians/150;
      }
      
      customParams.orientation=  'side_step';

    // }else if( (dis_agent_to_walls < 2.5) && clearance - sceneEntities[i].radius <= 0){
    }else if( (dis_btwn_agents < (cross_width+2) ) && clearance - sceneEntities[i].radius <= 0){
      // sceneEntities[i].agent.rotation.z = sceneEntities[i].agent.rotation.z + 1.5708;

    }else{
        customParams.orientation=  'front';
      }
*/

  orientationConstraint(sceneEntities[i], clearance, dis_btwn_agents);

  i += 1;
}

//=============== END ==  orientation constraint ================ For one agent passing through walls scenario ==================================



 /*
  if(customParams.scenario != 'bottleneck')
  {
//=========================================== our Long-range call started ======================================================================
    //Capsule to Capsule long-range collision avoidance.
    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {
        let [bestA, bestB, agent_i, agent_j] = getBestPoint(sceneEntities[i].x, sceneEntities[i].z, sceneEntities[j].x, sceneEntities[j].z);
        let [p_bestA, p_bestB, p_agent_i,p_agent_j] = getBestPoint(sceneEntities[i].px, sceneEntities[i].pz, sceneEntities[j].px, sceneEntities[j].pz);

        let [delta_correction_i, delta_correction_j, grad_i, grad_j, s] = longRangeConstraintCapsule(
            bestA, bestB,
            p_bestA, p_bestB,
            sceneEntities[i].agent.rotation.z, sceneEntities[j].agent.rotation.z,
            agent_i, agent_j,
            sceneEntities[i], sceneEntities[j],
            i, j
        );

        let long_stif = 1;

        sceneEntities[i].px += delta_correction_i.x;
        sceneEntities[i].pz += delta_correction_i.y;
        sceneEntities[j].px += delta_correction_j.x;
        sceneEntities[j].pz += delta_correction_j.y;

        // for utilities
        sceneEntities[i].grad.x += grad_i[0];
        sceneEntities[i].grad.z += grad_i[1];
        sceneEntities[j].grad.x += grad_j[0];
        sceneEntities[j].grad.z += grad_j[1];

        sceneEntities[i].grad.s = s;
        sceneEntities[j].grad.s = s;

        sceneEntities[i].grad.dx += long_stif * delta_correction_i.x;
        sceneEntities[i].grad.dz += long_stif * delta_correction_i.y;
        sceneEntities[j].grad.dx += long_stif * delta_correction_j.x;
        sceneEntities[j].grad.dz += long_stif * delta_correction_j.y;

        customParams.best[i][j] = [bestA, bestB]
        customParams.best[j][i] = [bestB, bestA]

        j += 1;
      }
      i += 1;
    }
//======================================== our Long-range call ended ===============================================================
  }
*/
  

//short range constraint.
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


    //Rotation constraint
    if(customParams.orientation === 'front')
    {
      i = 0;
      while (i < sceneEntities.length) {
      
        if(customParams.tempcount > 2) //allowing to adjust capsule's facing direction to the goal direction at the beginning.
        {
          let corrected_position = rotationConstraint(sceneEntities[i]);
        
          sceneEntities[i].px = corrected_position.x;
          sceneEntities[i].pz = corrected_position.z;
        }
      
        i += 1;
      }
  }

    pbdIters += 1;
  }
 

  sceneEntities.forEach(function (item) {

    const dx = item.px - item.x;
    const dz = item.pz - item.z;


    let cur_orientation = item.agent.rotation.z;

    if(customParams.orientation === 'front')
    {
        // item.agent.rotation.z = Math.atan2(dz, dx);  //original code for rotation
      let  next_orientation = Math.atan2(dz, dx);

        if( Math.abs(cur_orientation - next_orientation) >= 0.1)   // 0.0 makes the orientation better after recovering from the full body rotations.
      {
        item.agent.rotation.z = cur_orientation - cur_orientation/200;
      }
    }


// if(customParams.orientation === 'front')
// {
//   item.agent.rotation.z = Math.atan2(dz, dx);  //original code for rotation
// }
 
    item.vx = (item.px - item.x) / timestep;
    item.vz = (item.pz - item.z) / timestep;
    item.vy = (item.py - item.y) / timestep;

    item.x = item.px;
    item.z = item.pz;
    item.y = item.py;

  });

}
