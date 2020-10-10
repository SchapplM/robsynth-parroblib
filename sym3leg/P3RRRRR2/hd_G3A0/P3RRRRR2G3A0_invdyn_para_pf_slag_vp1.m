% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:38
% EndTime: 2020-03-09 21:11:44
% DurationCPUTime: 6.59s
% Computational Cost: add. (13233->499), mult. (32415->924), div. (10368->11), fcn. (30276->60), ass. (0->373)
t256 = xDP(2);
t246 = cos(qJ(3,1));
t215 = 0.1e1 / t246 ^ 2;
t238 = sin(qJ(2,1));
t206 = 0.1e1 / t238;
t262 = 0.1e1 / pkin(2);
t264 = 0.1e1 / pkin(1);
t336 = t262 * t264;
t292 = t206 * t336;
t273 = t215 * t292;
t227 = legFrame(1,2);
t192 = sin(t227);
t195 = cos(t227);
t239 = sin(qJ(1,1));
t247 = cos(qJ(2,1));
t248 = cos(qJ(1,1));
t138 = t238 * t239 - t248 * t247;
t213 = t246 ^ 2;
t326 = pkin(2) * t138 * t213;
t237 = sin(qJ(3,1));
t414 = pkin(1) * t247;
t329 = t237 * t414;
t353 = t195 * t237;
t413 = pkin(1) * t248;
t78 = -t192 * t326 + (-pkin(2) * t353 + t192 * t413) * t246 - t195 * t329;
t72 = t78 * t256 * t273;
t257 = xDP(1);
t356 = t192 * t237;
t81 = t195 * t326 + (-pkin(2) * t356 - t195 * t413) * t246 - t192 * t329;
t75 = t81 * t257 * t273;
t111 = pkin(2) * (t238 * t248 + t239 * t247) * t246 + t239 * pkin(1);
t255 = xDP(3);
t214 = 0.1e1 / t246;
t274 = t214 * t292;
t93 = t111 * t255 * t274;
t54 = t75 + t72 + t93;
t374 = t138 * t246;
t101 = t192 * t374 + t353;
t102 = -t195 * t374 + t356;
t224 = qJ(1,1) + qJ(2,1);
t174 = sin(t224);
t350 = t206 * t264;
t293 = t214 * t350;
t337 = t255 * t264;
t66 = -t174 * t206 * t337 + (t101 * t256 + t102 * t257) * t293;
t407 = -t54 - t66;
t442 = t407 ^ 2;
t243 = cos(qJ(3,2));
t212 = 0.1e1 / t243 ^ 2;
t235 = sin(qJ(2,2));
t205 = 0.1e1 / t235;
t294 = t205 * t336;
t275 = t212 * t294;
t226 = legFrame(2,2);
t191 = sin(t226);
t194 = cos(t226);
t236 = sin(qJ(1,2));
t244 = cos(qJ(2,2));
t245 = cos(qJ(1,2));
t137 = t235 * t236 - t245 * t244;
t210 = t243 ^ 2;
t327 = pkin(2) * t137 * t210;
t234 = sin(qJ(3,2));
t416 = pkin(1) * t244;
t330 = t234 * t416;
t354 = t194 * t234;
t415 = pkin(1) * t245;
t77 = -t191 * t327 + (-pkin(2) * t354 + t191 * t415) * t243 - t194 * t330;
t71 = t77 * t256 * t275;
t357 = t191 * t234;
t80 = t194 * t327 + (-pkin(2) * t357 - t194 * t415) * t243 - t191 * t330;
t74 = t80 * t257 * t275;
t110 = pkin(2) * (t235 * t245 + t236 * t244) * t243 + t236 * pkin(1);
t211 = 0.1e1 / t243;
t276 = t211 * t294;
t92 = t110 * t255 * t276;
t53 = t74 + t71 + t92;
t375 = t137 * t243;
t100 = -t194 * t375 + t357;
t221 = qJ(1,2) + qJ(2,2);
t171 = sin(t221);
t351 = t205 * t264;
t295 = t211 * t351;
t99 = t191 * t375 + t354;
t65 = -t171 * t205 * t337 + (t100 * t257 + t256 * t99) * t295;
t408 = -t53 - t65;
t441 = t408 ^ 2;
t240 = cos(qJ(3,3));
t209 = 0.1e1 / t240 ^ 2;
t232 = sin(qJ(2,3));
t204 = 0.1e1 / t232;
t296 = t204 * t336;
t277 = t209 * t296;
t225 = legFrame(3,2);
t190 = sin(t225);
t193 = cos(t225);
t233 = sin(qJ(1,3));
t241 = cos(qJ(2,3));
t242 = cos(qJ(1,3));
t136 = t232 * t233 - t242 * t241;
t207 = t240 ^ 2;
t328 = pkin(2) * t136 * t207;
t231 = sin(qJ(3,3));
t418 = pkin(1) * t241;
t331 = t231 * t418;
t355 = t193 * t231;
t417 = pkin(1) * t242;
t76 = -t190 * t328 + (-pkin(2) * t355 + t190 * t417) * t240 - t193 * t331;
t70 = t76 * t256 * t277;
t358 = t190 * t231;
t79 = t193 * t328 + (-pkin(2) * t358 - t193 * t417) * t240 - t190 * t331;
t73 = t79 * t257 * t277;
t109 = pkin(2) * (t232 * t242 + t233 * t241) * t240 + t233 * pkin(1);
t208 = 0.1e1 / t240;
t278 = t208 * t296;
t91 = t109 * t255 * t278;
t52 = t73 + t70 + t91;
t218 = qJ(1,3) + qJ(2,3);
t168 = sin(t218);
t352 = t204 * t264;
t297 = t208 * t352;
t376 = t136 * t240;
t97 = t190 * t376 + t355;
t98 = -t193 * t376 + t358;
t64 = -t168 * t204 * t337 + (t256 * t97 + t257 * t98) * t297;
t409 = -t52 - t64;
t440 = t409 ^ 2;
t439 = 2 * m(3);
t426 = m(3) * rSges(3,2);
t165 = rSges(3,1) * t426 - Icges(3,4);
t260 = 0.2e1 * qJ(3,1);
t199 = sin(t260);
t202 = cos(t260);
t431 = rSges(3,2) ^ 2;
t432 = rSges(3,1) ^ 2;
t151 = (-t431 + t432) * m(3) - Icges(3,1) + Icges(3,2);
t419 = t151 / 0.2e1;
t438 = -t165 * t199 + t202 * t419;
t259 = 0.2e1 * qJ(3,2);
t198 = sin(t259);
t201 = cos(t259);
t437 = -t165 * t198 + t201 * t419;
t258 = 0.2e1 * qJ(3,3);
t197 = sin(t258);
t200 = cos(t258);
t436 = -t165 * t197 + t200 * t419;
t425 = m(3) * rSges(3,3);
t158 = m(2) * rSges(2,2) - t425;
t254 = m(2) * rSges(2,1);
t270 = rSges(3,1) * t246 - rSges(3,2) * t237;
t435 = -t158 * t238 - (-t270 * m(3) - t254) * t247;
t271 = rSges(3,1) * t243 - rSges(3,2) * t234;
t434 = -t158 * t235 - (-t271 * m(3) - t254) * t244;
t272 = rSges(3,1) * t240 - rSges(3,2) * t231;
t433 = -t158 * t232 - (-t272 * m(3) - t254) * t241;
t58 = t64 ^ 2;
t59 = t65 ^ 2;
t60 = t66 ^ 2;
t430 = -2 * t151;
t429 = m(3) * pkin(1);
t428 = (m(1) * rSges(1,2));
t427 = m(3) * rSges(3,1);
t424 = rSges(2,2) * g(3);
t423 = rSges(3,2) * g(3);
t422 = pkin(1) * t58;
t421 = pkin(1) * t59;
t420 = pkin(1) * t60;
t25 = t73 / 0.2e1 + t70 / 0.2e1 + t91 / 0.2e1 + t64;
t412 = t52 * t25;
t26 = t74 / 0.2e1 + t71 / 0.2e1 + t92 / 0.2e1 + t65;
t411 = t53 * t26;
t27 = t75 / 0.2e1 + t72 / 0.2e1 + t93 / 0.2e1 + t66;
t410 = t54 * t27;
t139 = g(1) * t193 - g(2) * t190;
t406 = rSges(3,1) * t139;
t140 = g(1) * t194 - g(2) * t191;
t405 = rSges(3,1) * t140;
t141 = g(1) * t195 - g(2) * t192;
t404 = rSges(3,1) * t141;
t349 = t208 * t262;
t106 = (t190 * t257 + t193 * t256) * t349;
t403 = t106 * t409;
t347 = t211 * t262;
t107 = (t191 * t257 + t194 * t256) * t347;
t402 = t107 * t408;
t345 = t214 * t262;
t108 = (t192 * t257 + t195 * t256) * t345;
t401 = t108 * t407;
t103 = t106 ^ 2;
t400 = (t103 / 0.2e1 + t412) * t232;
t104 = t107 ^ 2;
t399 = (t104 / 0.2e1 + t411) * t235;
t105 = t108 ^ 2;
t398 = (t105 / 0.2e1 + t410) * t238;
t397 = t207 * t409;
t163 = -rSges(3,2) * t425 + Icges(3,6);
t164 = rSges(3,1) * t425 - Icges(3,5);
t118 = t163 * t240 - t164 * t231;
t309 = t109 * t349;
t288 = pkin(1) * t232 + rSges(3,3);
t94 = (-t288 * t426 + Icges(3,6)) * t240 - t231 * (t288 * t427 - Icges(3,5));
t396 = (t118 * t309 - t168 * t94) * t297;
t395 = t208 * t97;
t394 = t208 * t98;
t393 = t210 * t408;
t119 = t163 * t243 - t164 * t234;
t308 = t110 * t347;
t287 = pkin(1) * t235 + rSges(3,3);
t95 = (-t287 * t426 + Icges(3,6)) * t243 - t234 * (t287 * t427 - Icges(3,5));
t392 = (t119 * t308 - t171 * t95) * t295;
t391 = t211 * t99;
t390 = t213 * t407;
t120 = t163 * t246 - t164 * t237;
t307 = t111 * t345;
t286 = pkin(1) * t238 + rSges(3,3);
t96 = (-t286 * t426 + Icges(3,6)) * t246 - t237 * (t286 * t427 - Icges(3,5));
t389 = (t120 * t307 - t174 * t96) * t293;
t388 = t100 * t211;
t387 = t101 * t214;
t386 = t102 * t214;
t385 = t103 * t208;
t384 = t104 * t211;
t383 = t105 * t214;
t382 = t106 * t241;
t381 = t107 * t244;
t380 = t108 * t247;
t373 = t151 * t197;
t372 = t151 * t198;
t371 = t151 * t199;
t367 = t164 * t240;
t366 = t164 * t243;
t365 = t164 * t246;
t361 = t165 * t200;
t360 = t165 * t201;
t359 = t165 * t202;
t348 = t209 * t262;
t346 = t212 * t262;
t344 = t215 * t262;
t343 = t231 * t232;
t342 = t234 * t235;
t341 = t237 * t238;
t340 = t240 * t241;
t339 = t243 * t244;
t338 = t246 * t247;
t223 = -qJ(3,1) + qJ(2,1);
t222 = qJ(3,1) + qJ(2,1);
t220 = qJ(2,2) - qJ(3,2);
t219 = qJ(2,2) + qJ(3,2);
t217 = -qJ(3,3) + qJ(2,3);
t216 = qJ(3,3) + qJ(2,3);
t162 = m(1) * rSges(1,1) + m(2) * pkin(1);
t335 = (t162 + t429) * g(3);
t334 = 0.2e1 * pkin(1);
t333 = g(3) * t428;
t332 = t429 / 0.2e1;
t325 = t409 * t382;
t324 = t408 * t381;
t323 = t407 * t380;
t322 = t76 * t348;
t321 = t79 * t348;
t313 = t431 + t432;
t266 = Icges(2,3) + ((rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2)) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(3,3) ^ 2 + t313) * m(3)) / 0.2e1;
t82 = t266 + t436;
t320 = t82 * t348;
t319 = t77 * t346;
t318 = t80 * t346;
t83 = t266 + t437;
t317 = t83 * t346;
t316 = t78 * t344;
t315 = t81 * t344;
t84 = t266 + t438;
t314 = t84 * t344;
t312 = t106 * t343;
t311 = t107 * t342;
t310 = t108 * t341;
t306 = t118 * t348;
t305 = t119 * t346;
t304 = t120 * t344;
t303 = t190 * t349;
t302 = t191 * t347;
t301 = t192 * t345;
t300 = t193 * t349;
t299 = t194 * t347;
t298 = t195 * t345;
t285 = -pkin(1) * t427 / 0.2e1;
t196 = g(3) * t425;
t281 = (m(2) * (rSges(2,1) * t139 - t424) + t196) * t168;
t280 = (m(2) * (rSges(2,1) * t140 - t424) + t196) * t171;
t279 = (m(2) * (rSges(2,1) * t141 - t424) + t196) * t174;
t269 = t158 * t241 + t232 * t254;
t268 = t158 * t244 + t235 * t254;
t267 = t158 * t247 + t238 * t254;
t265 = Icges(1,3) + (m(3) + m(2)) * pkin(1) ^ 2 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t266;
t261 = pkin(2) ^ 2;
t252 = rSges(2,1) * g(3);
t251 = rSges(3,1) * g(3);
t230 = xDDP(1);
t229 = xDDP(2);
t228 = xDDP(3);
t189 = qJ(1,1) + t223;
t188 = qJ(1,1) + t222;
t187 = qJ(1,2) + t220;
t186 = qJ(1,2) + t219;
t185 = qJ(1,3) + t217;
t184 = qJ(1,3) + t216;
t183 = cos(t224);
t182 = cos(t223);
t181 = cos(t222);
t180 = cos(t221);
t179 = cos(t220);
t178 = cos(t219);
t177 = cos(t218);
t176 = cos(t217);
t175 = cos(t216);
t173 = sin(t223);
t172 = sin(t222);
t170 = sin(t220);
t169 = sin(t219);
t167 = sin(t217);
t166 = sin(t216);
t156 = t313 * m(3) + Icges(3,3);
t129 = rSges(3,2) * t141;
t128 = rSges(3,2) * t140;
t127 = rSges(3,2) * t139;
t117 = m(2) * (rSges(2,2) * t141 + t252);
t116 = m(2) * (rSges(2,2) * t140 + t252);
t115 = m(2) * (rSges(2,2) * t139 + t252);
t69 = t435 * pkin(1) + t84;
t68 = t434 * pkin(1) + t83;
t67 = t433 * pkin(1) + t82;
t57 = t435 * t334 + t265 + t438;
t56 = t434 * t334 + t265 + t437;
t55 = t433 * t334 + t265 + t436;
t51 = t156 * t301 + (t81 * t304 + t96 * t386) * t350;
t50 = t156 * t298 + (t78 * t304 + t96 * t387) * t350;
t49 = t156 * t302 + (t80 * t305 + t95 * t388) * t351;
t48 = t156 * t299 + (t77 * t305 + t95 * t391) * t351;
t47 = t156 * t303 + (t79 * t306 + t94 * t394) * t352;
t46 = t156 * t300 + (t76 * t306 + t94 * t395) * t352;
t45 = (-t174 * t69 + t84 * t307) * t350;
t44 = (-t171 * t68 + t83 * t308) * t351;
t43 = (-t168 * t67 + t82 * t309) * t352;
t42 = (-t174 * t57 + t69 * t307) * t350;
t41 = (-t171 * t56 + t68 * t308) * t351;
t40 = (-t168 * t55 + t67 * t309) * t352;
t33 = t120 * t301 + (t81 * t314 + t69 * t386) * t350;
t32 = t120 * t298 + (t78 * t314 + t69 * t387) * t350;
t31 = t119 * t302 + (t80 * t317 + t68 * t388) * t351;
t30 = t119 * t299 + (t77 * t317 + t68 * t391) * t351;
t29 = t118 * t303 + (t79 * t320 + t67 * t394) * t352;
t28 = t118 * t300 + (t76 * t320 + t67 * t395) * t352;
t24 = t96 * t301 + (t69 * t315 + t57 * t386) * t350;
t23 = t96 * t298 + (t69 * t316 + t57 * t387) * t350;
t22 = t95 * t302 + (t68 * t318 + t56 * t388) * t351;
t21 = t95 * t299 + (t68 * t319 + t56 * t391) * t351;
t20 = t94 * t303 + (t67 * t321 + t55 * t394) * t352;
t19 = t94 * t300 + (t67 * t322 + t55 * t395) * t352;
t15 = (-t60 * t414 + (-t246 * t442 - t383) * pkin(2)) * t350;
t14 = (-t59 * t416 + (-t243 * t441 - t384) * pkin(2)) * t351;
t13 = (-t58 * t418 + (-t240 * t440 - t385) * pkin(2)) * t352;
t12 = (-t261 * t390 + (pkin(1) * t66 + (0.2e1 * t27 * t338 - t310) * pkin(2)) * pkin(1)) * t66 * t274 + (-pkin(2) * t390 + (-t338 * t407 - t310) * pkin(1)) * t54 * t293 + ((pkin(1) * t341 * t407 + pkin(2) * t108) * t246 + pkin(1) * t380) * t215 * t108 * t350;
t11 = (-t261 * t393 + (pkin(1) * t65 + (0.2e1 * t26 * t339 - t311) * pkin(2)) * pkin(1)) * t65 * t276 + (-pkin(2) * t393 + (-t339 * t408 - t311) * pkin(1)) * t53 * t295 + ((pkin(1) * t342 * t408 + pkin(2) * t107) * t243 + pkin(1) * t381) * t212 * t107 * t351;
t10 = (-t261 * t397 + (pkin(1) * t64 + (0.2e1 * t25 * t340 - t312) * pkin(2)) * pkin(1)) * t64 * t278 + (-pkin(2) * t397 + (-t340 * t409 - t312) * pkin(1)) * t52 * t297 + ((pkin(1) * t343 * t409 + pkin(2) * t106) * t240 + pkin(1) * t382) * t209 * t106 * t352;
t9 = -t96 * t15 - t120 * t12 + t156 * t237 * t383 + m(3) * (-(g(1) * t192 + g(2) * t195) * t270 + (-g(3) * t174 + t141 * t183) * (rSges(3,1) * t237 + rSges(3,2) * t246)) + (t371 / 0.2e1 + t359) * t442 + (t173 * t285 + (rSges(3,1) * t172 + (t181 + t182) * rSges(3,2)) * t332) * t60;
t8 = -t95 * t14 - t119 * t11 + t156 * t234 * t384 + m(3) * (-(g(1) * t191 + g(2) * t194) * t271 + (-g(3) * t171 + t140 * t180) * (rSges(3,1) * t234 + rSges(3,2) * t243)) + (t372 / 0.2e1 + t360) * t441 + (t170 * t285 + (rSges(3,1) * t169 + (t178 + t179) * rSges(3,2)) * t332) * t59;
t7 = -t94 * t13 - t118 * t10 + t156 * t231 * t385 + m(3) * (-(g(1) * t190 + g(2) * t193) * t272 + (-g(3) * t168 + t139 * t177) * (rSges(3,1) * t231 + rSges(3,2) * t240)) + (t373 / 0.2e1 + t361) * t440 + (t167 * t285 + (rSges(3,1) * t166 + (t175 + t176) * rSges(3,2)) * t332) * t58;
t6 = -t69 * t15 - t84 * t12 + t117 * t183 + t279 - (-0.2e1 * t359 - t371) * t401 + t267 * t420 + (-t365 + (t120 * t214 - t163) * t237) * t105 + (t270 * t183 * g(3) + (-rSges(3,3) * t183 + t270 * t174) * t141 + ((-t182 / 0.2e1 + t181 / 0.2e1) * rSges(3,2) + (t173 / 0.2e1 + t172 / 0.2e1) * rSges(3,1)) * t420) * m(3);
t5 = -t68 * t14 - t83 * t11 + t116 * t180 + t280 - (-0.2e1 * t360 - t372) * t402 + t268 * t421 + (-t366 + (t119 * t211 - t163) * t234) * t104 + (t271 * t180 * g(3) + (-rSges(3,3) * t180 + t271 * t171) * t140 + ((-t179 / 0.2e1 + t178 / 0.2e1) * rSges(3,2) + (t170 / 0.2e1 + t169 / 0.2e1) * rSges(3,1)) * t421) * m(3);
t4 = -t67 * t13 - t82 * t10 + t115 * t177 + t281 - (-0.2e1 * t361 - t373) * t403 + t269 * t422 + (-t367 + (t118 * t208 - t163) * t231) * t103 + (t272 * t177 * g(3) + (-rSges(3,3) * t177 + t272 * t168) * t139 + ((-t176 / 0.2e1 + t175 / 0.2e1) * rSges(3,2) + (t167 / 0.2e1 + t166 / 0.2e1) * rSges(3,1)) * t422) * m(3);
t3 = -t57 * t15 - t69 * t12 + t239 * (t141 * t162 - t333) - (t237 * t246 * t430 + (-0.4e1 * t213 + 0.2e1) * t165) * t401 + (-t365 + (t214 * t96 - t163) * t237) * t105 + (-(t129 - t251) * cos(t189) / 0.2e1 - (-t404 - t423) * sin(t189) / 0.2e1 + (t129 + t251) * cos(t188) / 0.2e1 + (t404 - t423) * sin(t188) / 0.2e1) * m(3) + (-t141 * t425 + t117) * t183 + t279 + (t141 * t428 + t335) * t248 + (t239 * m(3) * t141 - 0.2e1 * t267 * t410 + ((-rSges(3,1) * t398 + rSges(3,2) * t323) * t246 + (rSges(3,1) * t323 + rSges(3,2) * t398) * t237) * t439) * pkin(1);
t2 = -t56 * t14 - t68 * t11 + t236 * (t140 * t162 - t333) - (t234 * t243 * t430 + (-0.4e1 * t210 + 0.2e1) * t165) * t402 + (-t366 + (t211 * t95 - t163) * t234) * t104 + (-(t128 - t251) * cos(t187) / 0.2e1 - (-t405 - t423) * sin(t187) / 0.2e1 + (t128 + t251) * cos(t186) / 0.2e1 + (t405 - t423) * sin(t186) / 0.2e1) * m(3) + (-t140 * t425 + t116) * t180 + t280 + (t140 * t428 + t335) * t245 + (t236 * m(3) * t140 - 0.2e1 * t268 * t411 + ((-rSges(3,1) * t399 + rSges(3,2) * t324) * t243 + (rSges(3,1) * t324 + rSges(3,2) * t399) * t234) * t439) * pkin(1);
t1 = -t55 * t13 - t67 * t10 + t233 * (t139 * t162 - t333) - (t231 * t240 * t430 + (-0.4e1 * t207 + 0.2e1) * t165) * t403 + (-t367 + (t208 * t94 - t163) * t231) * t103 + (-(t127 - t251) * cos(t185) / 0.2e1 - (-t406 - t423) * sin(t185) / 0.2e1 + (t127 + t251) * cos(t184) / 0.2e1 + (t406 - t423) * sin(t184) / 0.2e1) * m(3) + (-t139 * t425 + t115) * t177 + t281 + (t139 * t428 + t335) * t242 + (t233 * m(3) * t139 - 0.2e1 * t269 * t412 + ((-rSges(3,1) * t400 + rSges(3,2) * t325) * t240 + (rSges(3,1) * t325 + rSges(3,2) * t400) * t231) * t439) * pkin(1);
t16 = [(-g(1) + t230) * m(4) + ((t195 * t229 * t51 + (t230 * t51 + t9) * t192) * t214 + (t194 * t229 * t49 + (t230 * t49 + t8) * t191) * t211 + (t193 * t229 * t47 + (t230 * t47 + t7) * t190) * t208) * t262 + (((t24 * t386 + t33 * t315) * t230 + (t24 * t387 + t33 * t316) * t229 + (-t174 * t24 + t33 * t307) * t228 + t3 * t386 + t6 * t315) * t206 + ((t22 * t388 + t31 * t318) * t230 + (t22 * t391 + t31 * t319) * t229 + (-t171 * t22 + t31 * t308) * t228 + t2 * t388 + t5 * t318) * t205 + ((t20 * t394 + t29 * t321) * t230 + (t20 * t395 + t29 * t322) * t229 + (-t168 * t20 + t29 * t309) * t228 + t1 * t394 + t4 * t321) * t204) * t264; (-g(2) + t229) * m(4) + ((t192 * t230 * t50 + (t229 * t50 + t9) * t195) * t214 + (t191 * t230 * t48 + (t229 * t48 + t8) * t194) * t211 + (t190 * t230 * t46 + (t229 * t46 + t7) * t193) * t208) * t262 + (((t23 * t386 + t32 * t315) * t230 + (t23 * t387 + t32 * t316) * t229 + (-t174 * t23 + t32 * t307) * t228 + t3 * t387 + t6 * t316) * t206 + ((t21 * t388 + t30 * t318) * t230 + (t21 * t391 + t30 * t319) * t229 + (-t171 * t21 + t30 * t308) * t228 + t2 * t391 + t5 * t319) * t205 + ((t19 * t394 + t28 * t321) * t230 + (t19 * t395 + t28 * t322) * t229 + (-t168 * t19 + t28 * t309) * t228 + t1 * t395 + t4 * t322) * t204) * t264; (-g(3) + t228) * m(4) + ((t190 * t396 + t191 * t392 + t192 * t389) * t230 + (t193 * t396 + t194 * t392 + t195 * t389) * t229) * t262 + (((t45 * t315 + t42 * t386) * t230 + (t45 * t316 + t42 * t387) * t229 + (-t174 * t42 + t45 * t307) * t228 - t174 * t3 + t6 * t307) * t206 + ((t44 * t318 + t41 * t388) * t230 + (t44 * t319 + t41 * t391) * t229 + (-t171 * t41 + t44 * t308) * t228 - t171 * t2 + t5 * t308) * t205 + ((t43 * t321 + t40 * t394) * t230 + (t43 * t322 + t40 * t395) * t229 + (-t168 * t40 + t43 * t309) * t228 - t168 * t1 + t4 * t309) * t204) * t264;];
tauX  = t16;
