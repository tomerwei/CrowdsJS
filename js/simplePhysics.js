export function step(RADIUS, sceneEntities, world, customParams = {}) {
  // from https://wickedengine.net/2020/04/26/capsule-collision-detection/
  function ClosestPointOnLineSegment(A, B, Point) {
    const AB = B - A;
    const t = AB.dot(Point - A) / AB.dot(AB);
    return A + Math.min((Math.max(t, 0), 1)) * AB;
  }

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

  function distance(x1, y1, x2, y2) {
    return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
  }

  const AGENTSIZE = RADIUS * 2;

  const epsilon = 0.0001;
  const timestep = 0.05;

  // collision functions
  function is_colliding(x11, y11, x12, y12, x21, y21, x22, y22) {
    // console.log(segments_distance(x11, y11, x12, y12, x21, y21, x22, y22));
    return (
      segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) < RADIUS * 2
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

  var i = 0,
    j;
  while (i < sceneEntities.length) {
    j = i + 1;
    while (j < sceneEntities.length) {
      /*  TODO modify lines below  */
      /*  -----------------------  */
      const iCoords = rotateLineSegment(
        sceneEntities[i].agent.position.x,
        sceneEntities[i].agent.position.z - agentLength,
        sceneEntities[i].agent.position.x,
        sceneEntities[i].agent.position.z + agentLength,
        sceneEntities[i].agent.rotation.z
      );
      const jCoords = rotateLineSegment(
        sceneEntities[j].agent.position.x,
        sceneEntities[j].agent.position.z - agentLength,
        sceneEntities[j].agent.position.x,
        sceneEntities[j].agent.position.z + agentLength,
        sceneEntities[j].agent.rotation.z
      );
      // if (
      //   is_colliding(
      //     iCoords[0],
      //     iCoords[1],
      //     iCoords[2],
      //     iCoords[3],
      //     jCoords[0],
      //     jCoords[1],
      //     jCoords[2],
      //     jCoords[3]
      //   )
      // ) {
      //   sceneEntities[i].colliding = true;
      //   sceneEntities[j].colliding = true;
      // }

      // collisions from pdf/website
      // Agent A
      const a = {
        tip: new THREE.Vector2(iCoords[0], iCoords[1]),
        base: new THREE.Vector2(iCoords[2], iCoords[3]),
        radius: RADIUS,
      };
      // Agent B
      const b = {
        tip: new THREE.Vector2(jCoords[0], jCoords[1]),
        base: new THREE.Vector2(jCoords[2], jCoords[3]),
        radius: RADIUS,
      };

      // capsule A:
      const a_Normal = a.tip.normalize() - a.base.normalize();
      const a_LineEndOffset = a_Normal * a.radius;
      const a_A = a.base + a_LineEndOffset;
      const a_B = a.tip - a_LineEndOffset;

      // capsule B:
      const b_Normal = b.tip.normalize() - b.base.normalize();
      const b_LineEndOffset = b_Normal * b.radius;
      const b_A = b.base + b_LineEndOffset;
      const b_B = b.tip + b_LineEndOffset;

      // vectors between line endpoints:
      const v0 = b_A - a_A;
      const v1 = b_B - a_A;
      const v2 = b_A - a_B;
      const v3 = b_B - a_B;

      // squared distances:
      const d0 = v0.dot(v0);
      const d1 = v1.dot(v1);
      const d2 = v2.dot(v2);
      const d3 = v3.dot(v3);

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

      const penetration_normal = bestA - bestB;
      const len = penetration_normal.length;
      penetration_normal /= len; // normalize
      const penetration_depth = a.radius + b.radius - len;
      const intersects = penetration_depth > 0;

      if (intersects) {
        sceneEntities[i].colliding = true;
        sceneEntities[j].colliding = true;

        // this needs to go below
        // sceneEntities[i].agent.position.x +=
        //   penetration_normal.x * 0.5 * penetration_normal.x;
        // sceneEntities[i].agent.position.z +=
        //   penetration_normal.z * 0.5 * penetration_normal.y;

        // sceneEntities[j].agent.position.x +=
        //   -1 * penetration_normal.x * 0.5 * penetration_normal.x;
        // sceneEntities[j].agent.position.z +=
        //   -1 * penetration_normal.z * 0.5 * penetration_normal.y;

        sceneEntities[i].cvx +=
          penetration_normal.x * 0.5 * penetration_normal.x;
        sceneEntities[i].cvz +=
          penetration_normal.z * 0.5 * penetration_normal.y;

        sceneEntities[j].cvx +=
          -1 * penetration_normal.x * 0.5 * penetration_normal.x;
        sceneEntities[j].cvz +=
          -1 * penetration_normal.z * 0.5 * penetration_normal.y;
      }

      // collision avoidance from https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
      /*
      const d = distance(
        sceneEntities[i].agent.position.x,
        sceneEntities[i].agent.position.z,
        sceneEntities[j].agent.position.x,
        sceneEntities[j].agent.position.z
      );
      const m1 = 10;
      const w1 = 1 / m1;
      const m2 = 10;
      const w2 = 1 / m2;

      const p1x = sceneEntities[i].agent.position.x;
      const p1y = sceneEntities[i].agent.position.z;
      const p2x = sceneEntities[j].agent.position.x;
      const p2y = sceneEntities[j].agent.position.z;

      const dx =
        (w1 / (w1 + w2)) *
        (Math.abs(p1x - p2x) - d) *
        ((p1x - p2x) / Math.abs(p1x - p2x));

      const dy =
        (w1 / (w1 + w2)) *
        (Math.abs(p1y - p2y) - d) *
        ((p1y - p2y) / Math.abs(p1y - p2y));

      sceneEntities[i].agent.position.x += -dx;
      sceneEntities[i].agent.position.z += -dy;

      sceneEntities[j].agent.position.x += dx;
      sceneEntities[j].agent.position.z += dy;
      */
      /*  -----------------------  */
      j += 1;
    }
    i += 1;
  }

  console.log("Length:", sceneEntities.length);
  sceneEntities.forEach(function (item) {
    // if (item.colliding) {
    //   item.x += timestep * item.cvx;
    //   item.z += timestep * item.cvz;
    // } else {
    //   item.x += timestep * item.vx;
    //   item.z += timestep * item.vz;
    // }

    if (item.x < -world.x / 2) {
      item.x = world.x / 2 - epsilon;
    } else if (item.x > world.x / 2) {
      item.x = -world.x / 2 + epsilon;
    }
    if (item.z < -world.z / 2) {
      item.z = world.z / 2 - epsilon;
    } else if (item.z > world.z / 2) {
      item.z = -world.z / 2 + epsilon;
    }
  });
}
