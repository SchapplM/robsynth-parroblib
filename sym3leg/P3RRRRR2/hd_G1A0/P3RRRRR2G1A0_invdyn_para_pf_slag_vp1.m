% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR2G1P1A0
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:01
% EndTime: 2020-03-09 21:04:07
% DurationCPUTime: 6.24s
% Computational Cost: add. (16938->475), mult. (17022->829), div. (7512->18), fcn. (11136->84), ass. (0->384)
t458 = 2 * pkin(1);
t457 = -2 * rSges(3,1);
t456 = 2 * rSges(3,2);
t435 = m(3) * rSges(3,2);
t161 = rSges(3,1) * t435 - Icges(3,4);
t253 = 2 * qJ(3,1);
t197 = sin(t253);
t200 = cos(t253);
t443 = rSges(3,2) ^ 2;
t444 = rSges(3,1) ^ 2;
t136 = (-t443 + t444) * m(3) - Icges(3,1) + Icges(3,2);
t430 = t136 / 0.2e1;
t455 = -t161 * t197 + t200 * t430;
t252 = 2 * qJ(3,2);
t196 = sin(t252);
t199 = cos(t252);
t454 = -t161 * t196 + t199 * t430;
t251 = 2 * qJ(3,3);
t195 = sin(t251);
t198 = cos(t251);
t453 = -t161 * t195 + t198 * t430;
t434 = m(3) * rSges(3,3);
t157 = m(2) * rSges(2,2) - t434;
t237 = sin(qJ(2,1));
t243 = cos(qJ(2,1));
t247 = m(2) * rSges(2,1);
t236 = sin(qJ(3,1));
t242 = cos(qJ(3,1));
t449 = -rSges(3,1) * t242 + rSges(3,2) * t236;
t452 = -t157 * t237 - (m(3) * t449 - t247) * t243;
t235 = sin(qJ(2,2));
t241 = cos(qJ(2,2));
t234 = sin(qJ(3,2));
t240 = cos(qJ(3,2));
t448 = -rSges(3,1) * t240 + rSges(3,2) * t234;
t451 = -t157 * t235 - (m(3) * t448 - t247) * t241;
t233 = sin(qJ(2,3));
t239 = cos(qJ(2,3));
t232 = sin(qJ(3,3));
t238 = cos(qJ(3,3));
t447 = -rSges(3,1) * t238 + rSges(3,2) * t232;
t450 = -t157 * t233 - (m(3) * t447 - t247) * t239;
t446 = -2 * pkin(1);
t442 = -4 * t161;
t441 = 2 * t161;
t440 = -m(3) / 0.2e1;
t439 = m(3) / 0.2e1;
t438 = m(3) * pkin(1);
t437 = m(1) * rSges(1,2);
t436 = m(3) * rSges(3,1);
t248 = xDP(3);
t249 = xDP(2);
t250 = xDP(1);
t205 = 0.1e1 / t238;
t201 = 0.1e1 / t233;
t258 = 1 / pkin(1);
t367 = t201 * t258;
t312 = t232 * t367;
t281 = t205 * t312;
t226 = legFrame(3,3);
t162 = qJ(1,3) + t226;
t154 = qJ(2,3) + t162;
t145 = cos(t154);
t315 = t145 * t367;
t142 = sin(t154);
t318 = t142 * t367;
t67 = t248 * t281 + t249 * t318 + t250 * t315;
t64 = t67 ^ 2;
t433 = pkin(1) * t64;
t209 = 0.1e1 / t240;
t202 = 0.1e1 / t235;
t365 = t202 * t258;
t311 = t234 * t365;
t280 = t209 * t311;
t222 = qJ(2,2) + qJ(1,2);
t227 = legFrame(2,3);
t155 = t227 + t222;
t146 = cos(t155);
t314 = t146 * t365;
t143 = sin(t155);
t317 = t143 * t365;
t68 = t248 * t280 + t249 * t317 + t250 * t314;
t65 = t68 ^ 2;
t432 = pkin(1) * t65;
t213 = 0.1e1 / t242;
t203 = 0.1e1 / t237;
t363 = t203 * t258;
t310 = t236 * t363;
t279 = t213 * t310;
t228 = legFrame(1,3);
t164 = qJ(1,1) + t228;
t156 = qJ(2,1) + t164;
t147 = cos(t156);
t313 = t147 * t363;
t144 = sin(t156);
t316 = t144 * t363;
t69 = t248 * t279 + t249 * t316 + t250 * t313;
t66 = t69 ^ 2;
t431 = pkin(1) * t66;
t429 = pkin(1) * t233;
t428 = pkin(1) * t235;
t427 = pkin(1) * t237;
t426 = pkin(1) * t239;
t425 = pkin(1) * t241;
t424 = pkin(1) * t243;
t204 = t238 ^ 2;
t423 = pkin(2) * t204;
t208 = t240 ^ 2;
t422 = pkin(2) * t208;
t212 = t242 ^ 2;
t421 = pkin(2) * t212;
t420 = pkin(2) * t238;
t419 = pkin(2) * t240;
t418 = pkin(2) * t242;
t206 = 0.1e1 / t238 ^ 2;
t255 = 0.1e1 / pkin(2);
t361 = t206 * t255;
t321 = (t420 + t426) * t361;
t368 = t201 * t232;
t275 = t321 * t368;
t269 = t258 * t275;
t266 = t248 * t269;
t217 = qJ(2,3) + qJ(3,3);
t165 = sin(t217);
t218 = qJ(2,3) - qJ(3,3);
t166 = sin(t218);
t122 = 0.1e1 / (t165 + t166);
t347 = t255 * t258;
t324 = t122 * t347;
t148 = qJ(3,3) + t154;
t149 = -qJ(3,3) + t154;
t79 = sin(t162) * t446 + (-sin(t149) - sin(t148)) * pkin(2);
t287 = t79 * t324;
t58 = t249 * t287;
t82 = cos(t162) * t446 + (-cos(t149) - cos(t148)) * pkin(2);
t286 = t82 * t324;
t61 = t250 * t286;
t37 = t61 / 0.2e1 + t58 / 0.2e1 - t266 / 0.2e1 + t67;
t46 = t58 + t61 - t266;
t417 = t46 * t37;
t210 = 0.1e1 / t240 ^ 2;
t359 = t210 * t255;
t320 = (t419 + t425) * t359;
t366 = t202 * t234;
t274 = t320 * t366;
t268 = t258 * t274;
t265 = t248 * t268;
t220 = qJ(2,2) + qJ(3,2);
t168 = sin(t220);
t221 = qJ(2,2) - qJ(3,2);
t169 = sin(t221);
t123 = 0.1e1 / (t168 + t169);
t323 = t123 * t347;
t191 = qJ(1,2) + t220;
t150 = t227 + t191;
t192 = qJ(1,2) + t221;
t151 = t227 + t192;
t163 = qJ(1,2) + t227;
t80 = sin(t163) * t446 + (-sin(t151) - sin(t150)) * pkin(2);
t285 = t80 * t323;
t59 = t249 * t285;
t83 = cos(t163) * t446 + (-cos(t151) - cos(t150)) * pkin(2);
t284 = t83 * t323;
t62 = t250 * t284;
t38 = t62 / 0.2e1 + t59 / 0.2e1 - t265 / 0.2e1 + t68;
t47 = t59 + t62 - t265;
t416 = t47 * t38;
t214 = 0.1e1 / t242 ^ 2;
t357 = t214 * t255;
t319 = (t418 + t424) * t357;
t364 = t203 * t236;
t273 = t319 * t364;
t267 = t258 * t273;
t264 = t248 * t267;
t223 = qJ(2,1) + qJ(3,1);
t171 = sin(t223);
t224 = qJ(2,1) - qJ(3,1);
t172 = sin(t224);
t124 = 0.1e1 / (t171 + t172);
t322 = t124 * t347;
t152 = qJ(3,1) + t156;
t153 = -qJ(3,1) + t156;
t81 = sin(t164) * t446 + (-sin(t153) - sin(t152)) * pkin(2);
t283 = t81 * t322;
t60 = t249 * t283;
t84 = cos(t164) * t446 + (-cos(t153) - cos(t152)) * pkin(2);
t282 = t84 * t322;
t63 = t250 * t282;
t39 = t63 / 0.2e1 + t60 / 0.2e1 - t264 / 0.2e1 + t69;
t48 = t60 + t63 - t264;
t415 = t48 * t39;
t183 = sin(t226);
t186 = cos(t226);
t112 = -g(1) * t183 + g(2) * t186;
t414 = rSges(3,1) * t112;
t184 = sin(t227);
t187 = cos(t227);
t113 = -g(1) * t184 + g(2) * t187;
t413 = rSges(3,1) * t113;
t185 = sin(t228);
t188 = cos(t228);
t114 = -g(1) * t185 + g(2) * t188;
t412 = rSges(3,1) * t114;
t115 = g(1) * t186 + g(2) * t183;
t411 = rSges(3,1) * t115;
t116 = g(1) * t187 + g(2) * t184;
t410 = rSges(3,1) * t116;
t117 = g(1) * t188 + g(2) * t185;
t409 = rSges(3,1) * t117;
t402 = rSges(3,3) * t112;
t401 = rSges(3,3) * t113;
t400 = rSges(3,3) * t114;
t399 = rSges(3,3) * t115;
t398 = rSges(3,3) * t116;
t397 = rSges(3,3) * t117;
t393 = t122 * t255;
t392 = t123 * t255;
t391 = t124 * t255;
t390 = t136 * t195;
t389 = t136 * t196;
t388 = t136 * t197;
t386 = t142 * t201;
t385 = t143 * t202;
t384 = t144 * t203;
t383 = t145 * t201;
t382 = t146 * t202;
t381 = t147 * t203;
t159 = -rSges(3,2) * t434 + Icges(3,6);
t377 = t159 * t206;
t376 = t159 * t210;
t375 = t159 * t214;
t371 = t161 * t198;
t370 = t161 * t199;
t369 = t161 * t200;
t362 = t205 * t255;
t360 = t209 * t255;
t358 = t213 * t255;
t216 = t248 ^ 2;
t356 = t216 * t255;
t355 = t216 / pkin(2) ^ 2;
t229 = xDDP(3);
t354 = t229 * t258;
t230 = xDDP(2);
t353 = t230 * t258;
t231 = xDDP(1);
t352 = t231 * t258;
t351 = t232 * t233;
t350 = t234 * t235;
t349 = t236 * t237;
t348 = t248 * t255;
t345 = pkin(1) * t436;
t344 = pkin(1) * t435;
t343 = t438 / 0.2e1;
t342 = pkin(1) * t348;
t341 = t239 * t423;
t340 = t241 * t422;
t339 = t243 * t421;
t338 = t79 * t393;
t337 = t82 * t393;
t160 = rSges(3,1) * t434 - Icges(3,5);
t85 = t159 * t238 - t160 * t232;
t336 = t85 * t393;
t335 = t80 * t392;
t334 = t83 * t392;
t86 = t159 * t240 - t160 * t234;
t333 = t86 * t392;
t332 = t81 * t391;
t331 = t84 * t391;
t87 = t159 * t242 - t160 * t236;
t330 = t87 * t391;
t43 = t46 + t67;
t329 = t43 * t348;
t44 = t47 + t68;
t328 = t44 * t348;
t45 = t48 + t69;
t327 = t45 * t348;
t326 = -2 * t345;
t325 = t443 + t444;
t309 = t232 * t355;
t308 = t234 * t355;
t307 = t236 * t355;
t306 = t248 * t351;
t305 = t248 * t350;
t304 = t248 * t349;
t303 = t160 * t355;
t302 = m(1) * rSges(1,1) + m(2) * pkin(1) + t438;
t298 = t355 / 0.2e1;
t297 = rSges(3,3) + t429;
t296 = rSges(3,3) + t428;
t295 = rSges(3,3) + t427;
t294 = m(3) * (t206 * t298 + t417) * t429;
t293 = m(3) * (t210 * t298 + t416) * t428;
t292 = m(3) * (t214 * t298 + t415) * t427;
t291 = -t345 / 0.2e1;
t272 = t157 * t239 + t233 * t247;
t271 = t157 * t241 + t235 * t247;
t270 = t157 * t243 + t237 * t247;
t263 = Icges(2,3) + ((rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2)) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (2 * rSges(3,3) ^ 2 + t325) * t439;
t55 = t263 + t453;
t56 = t263 + t454;
t57 = t263 + t455;
t257 = pkin(1) ^ 2;
t262 = Icges(1,3) + ((m(3) + m(2)) * t257) + ((rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1)) + t263;
t261 = t242 * t212;
t260 = t240 * t208;
t259 = t238 * t204;
t254 = pkin(2) ^ 2;
t225 = qJ(1,1) + qJ(2,1);
t219 = qJ(1,3) + qJ(2,3);
t215 = 0.1e1 / t261;
t211 = 0.1e1 / t260;
t207 = 0.1e1 / t259;
t194 = qJ(1,1) + t224;
t193 = qJ(1,1) + t223;
t190 = qJ(1,3) + t218;
t189 = qJ(1,3) + t217;
t182 = cos(t225);
t181 = cos(t224);
t180 = cos(t223);
t179 = cos(t222);
t178 = cos(t221);
t177 = cos(t220);
t176 = cos(t219);
t175 = cos(t218);
t174 = cos(t217);
t173 = sin(t225);
t170 = sin(t222);
t167 = sin(t219);
t141 = t325 * m(3) + Icges(3,3);
t105 = rSges(3,2) * t117;
t104 = rSges(3,2) * t116;
t103 = rSges(3,2) * t115;
t102 = rSges(3,2) * t114;
t101 = rSges(3,2) * t113;
t100 = rSges(3,2) * t112;
t78 = (-t295 * t435 + Icges(3,6)) * t242 - t236 * (t295 * t436 - Icges(3,5));
t77 = (-t296 * t435 + Icges(3,6)) * t240 - t234 * (t296 * t436 - Icges(3,5));
t76 = (-t297 * t435 + Icges(3,6)) * t238 - t232 * (t297 * t436 - Icges(3,5));
t75 = m(2) * (rSges(2,1) * t117 + rSges(2,2) * t114);
t74 = m(2) * (rSges(2,1) * t116 + rSges(2,2) * t113);
t73 = m(2) * (rSges(2,1) * t115 + rSges(2,2) * t112);
t72 = m(2) * (-rSges(2,1) * t114 + rSges(2,2) * t117);
t71 = m(2) * (-rSges(2,1) * t113 + rSges(2,2) * t116);
t70 = m(2) * (-rSges(2,1) * t112 + rSges(2,2) * t115);
t54 = t452 * pkin(1) + t57;
t53 = t451 * pkin(1) + t56;
t52 = t450 * pkin(1) + t55;
t51 = t452 * t458 + t262 + t455;
t50 = t451 * t458 + t262 + t454;
t49 = t450 * t458 + t262 + t453;
t36 = t254 * t45 * t261;
t35 = t254 * t44 * t260;
t34 = t254 * t43 * t259;
t33 = t87 * t358 + (t213 * t54 - t57 * t319) * t310;
t32 = t86 * t360 + (t209 * t53 - t56 * t320) * t311;
t31 = t85 * t362 + (t205 * t52 - t55 * t321) * t312;
t30 = (t57 * t331 + t54 * t381) * t258;
t29 = (t56 * t334 + t53 * t382) * t258;
t28 = (t55 * t337 + t52 * t383) * t258;
t27 = (t57 * t332 + t54 * t384) * t258;
t26 = (t56 * t335 + t53 * t385) * t258;
t25 = (t55 * t338 + t52 * t386) * t258;
t24 = (t54 * t331 + t51 * t381) * t258;
t23 = (t53 * t334 + t50 * t382) * t258;
t22 = (t52 * t337 + t49 * t383) * t258;
t21 = (t54 * t332 + t51 * t384) * t258;
t20 = (t53 * t335 + t50 * t385) * t258;
t19 = (t52 * t338 + t49 * t386) * t258;
t18 = t78 * t358 + (t213 * t51 - t54 * t319) * t310;
t17 = t77 * t360 + (t209 * t50 - t53 * t320) * t311;
t16 = t76 * t362 + (t205 * t49 - t52 * t321) * t312;
t12 = ((-t242 * t69 * t424 - t45 * t421) * t213 * t69 - t45 * t48 * t418 - t215 * t356) * t363;
t11 = ((-t240 * t68 * t425 - t44 * t422) * t209 * t68 - t44 * t47 * t419 - t211 * t356) * t365;
t10 = ((-t238 * t67 * t426 - t43 * t423) * t205 * t67 - t43 * t46 * t420 - t207 * t356) * t367;
t9 = (((-pkin(1) * t45 * t349 + t213 * t248) * t242 + t243 * t213 * t342) * t215 * t348 + ((t36 + t39 * t339 * t458 + (-pkin(1) * t213 * t304 + t257 * t69) * t242) * t69 + (t36 + (t45 * t339 - t304) * pkin(1)) * t48) * t357) * t363;
t8 = (((-pkin(1) * t44 * t350 + t209 * t248) * t240 + t241 * t209 * t342) * t211 * t348 + ((t35 + t38 * t340 * t458 + (-pkin(1) * t209 * t305 + t257 * t68) * t240) * t68 + (t35 + (t44 * t340 - t305) * pkin(1)) * t47) * t359) * t365;
t7 = (((-pkin(1) * t43 * t351 + t205 * t248) * t238 + t239 * t205 * t342) * t207 * t348 + ((t34 + t37 * t341 * t458 + (-pkin(1) * t205 * t306 + t257 * t67) * t238) * t67 + (t34 + (t43 * t341 - t306) * pkin(1)) * t46) * t361) * t367;
t6 = -t54 * t12 + t75 * t173 + t72 * t182 - t57 * t9 + (t215 * t87 - t375) * t307 + t270 * t431 + (-t303 + (-0.2e1 * t369 - t388) * t327) * t213 + ((t449 * t114 - t397) * t182 + (-t449 * t117 - t400) * t173 + ((-t181 / 0.2e1 + t180 / 0.2e1) * rSges(3,2) + (t172 / 0.2e1 + t171 / 0.2e1) * rSges(3,1)) * t431) * m(3);
t5 = -t53 * t11 + t74 * t170 + t71 * t179 - t56 * t8 + (t211 * t86 - t376) * t308 + t271 * t432 + (-t303 + (-0.2e1 * t370 - t389) * t328) * t209 + ((t448 * t113 - t398) * t179 + (-t448 * t116 - t401) * t170 + ((-t178 / 0.2e1 + t177 / 0.2e1) * rSges(3,2) + (t169 / 0.2e1 + t168 / 0.2e1) * rSges(3,1)) * t432) * m(3);
t4 = -t52 * t10 + t73 * t167 + t70 * t176 - t55 * t7 + (t207 * t85 - t377) * t309 + t272 * t433 + (-t303 + (-0.2e1 * t371 - t390) * t329) * t205 + ((t447 * t112 - t399) * t176 + (-t447 * t115 - t402) * t167 + ((-t175 / 0.2e1 + t174 / 0.2e1) * rSges(3,2) + (t166 / 0.2e1 + t165 / 0.2e1) * rSges(3,1)) * t433) * m(3);
t3 = -t51 * t12 - t54 * t9 + (-t214 * t303 + t292 * t457) * t242 + sin(qJ(1,1)) * (t114 * t437 + t302 * t117) - t270 * t458 * t415 + (t292 * t456 + (t215 * t78 - t375) * t355) * t236 + (t242 * t442 + (0.2e1 * (-t136 * t236 - t243 * t344) * t242 + t243 * t236 * t326 + t441) * t213) * t327 + (-m(3) * t397 + t72) * t182 + (-m(3) * t400 + t75) * t173 + (-t302 * t114 + t117 * t437) * cos(qJ(1,1)) + ((t105 + t412) * cos(t194) + (t102 - t409) * sin(t194)) * t440 + ((t105 - t412) * cos(t193) + (t102 + t409) * sin(t193)) * t439;
t2 = -t50 * t11 - t53 * t8 + (-t210 * t303 + t293 * t457) * t240 + sin(qJ(1,2)) * (t113 * t437 + t302 * t116) - t271 * t458 * t416 + (t293 * t456 + (t211 * t77 - t376) * t355) * t234 + (t240 * t442 + (0.2e1 * (-t136 * t234 - t241 * t344) * t240 + t241 * t234 * t326 + t441) * t209) * t328 + (-m(3) * t398 + t71) * t179 + (-m(3) * t401 + t74) * t170 + (-t302 * t113 + t116 * t437) * cos(qJ(1,2)) + ((t104 + t413) * cos(t192) + (t101 - t410) * sin(t192)) * t440 + ((t104 - t413) * cos(t191) + (t101 + t410) * sin(t191)) * t439;
t1 = -t49 * t10 - t52 * t7 + (-t206 * t303 + t294 * t457) * t238 + sin(qJ(1,3)) * (t112 * t437 + t302 * t115) - t272 * t458 * t417 + (t294 * t456 + (t207 * t76 - t377) * t355) * t232 + (t238 * t442 + (0.2e1 * (-t136 * t232 - t239 * t344) * t238 + t239 * t232 * t326 + t441) * t205) * t329 + (-m(3) * t399 + t70) * t176 + (-m(3) * t402 + t73) * t167 + (-t302 * t112 + t115 * t437) * cos(qJ(1,3)) + ((t103 + t414) * cos(t190) + (t100 - t411) * sin(t190)) * t440 + ((t103 - t414) * cos(t189) + (t100 + t411) * sin(t189)) * t439;
t13 = [t1 * t315 + t2 * t314 + t6 * t282 + t5 * t284 + t4 * t286 + t3 * t313 + (-g(1) + t231) * m(4) + (-t30 * t273 + (t24 * t364 + (t84 * t330 + t78 * t381) * t255) * t213 - t29 * t274 + (t23 * t366 + (t83 * t333 + t77 * t382) * t255) * t209 - t28 * t275 + (t22 * t368 + (t82 * t336 + t76 * t383) * t255) * t205) * t354 + (t22 * t386 + t23 * t385 + t24 * t384 + t28 * t338 + t29 * t335 + t30 * t332) * t353 + (t22 * t383 + t23 * t382 + t24 * t381 + t28 * t337 + t29 * t334 + t30 * t331) * t352; t1 * t318 + t2 * t317 + t6 * t283 + t5 * t285 + t4 * t287 + t3 * t316 + (-g(2) + t230) * m(4) + (-t27 * t273 + (t21 * t364 + (t81 * t330 + t78 * t384) * t255) * t213 - t26 * t274 + (t20 * t366 + (t80 * t333 + t77 * t385) * t255) * t209 - t25 * t275 + (t19 * t368 + (t79 * t336 + t76 * t386) * t255) * t205) * t354 + (t19 * t386 + t20 * t385 + t21 * t384 + t25 * t338 + t26 * t335 + t27 * t332) * t353 + (t19 * t383 + t20 * t382 + t21 * t381 + t25 * t337 + t26 * t334 + t27 * t331) * t352; t3 * t279 - t6 * t267 + (-t78 * t12 - t87 * t9 + t141 * t215 * t307 + m(3) * (g(3) * t449 + (t114 * t173 + t117 * t182) * (rSges(3,1) * t236 + rSges(3,2) * t242)) + (t388 / 0.2e1 + t369) * t45 ^ 2 + (t172 * t291 + (t171 * rSges(3,1) + (t180 + t181) * rSges(3,2)) * t343) * t66) * t358 + t2 * t280 - t5 * t268 + (-t77 * t11 - t86 * t8 + t141 * t211 * t308 + m(3) * (g(3) * t448 + (t113 * t170 + t116 * t179) * (rSges(3,1) * t234 + rSges(3,2) * t240)) + (t389 / 0.2e1 + t370) * t44 ^ 2 + (t169 * t291 + (t168 * rSges(3,1) + (t177 + t178) * rSges(3,2)) * t343) * t65) * t360 + t1 * t281 - t4 * t269 + (-t76 * t10 - t85 * t7 + t141 * t207 * t309 + m(3) * (g(3) * t447 + (t112 * t167 + t115 * t176) * (rSges(3,1) * t232 + rSges(3,2) * t238)) + (t390 / 0.2e1 + t371) * t43 ^ 2 + (t166 * t291 + (t165 * rSges(3,1) + (t174 + t175) * rSges(3,2)) * t343) * t64) * t362 - m(4) * g(3) + (t16 * t386 + t17 * t385 + t18 * t384 + t31 * t338 + t32 * t335 + t33 * t332) * t353 + (t16 * t383 + t17 * t382 + t18 * t381 + t31 * t337 + t32 * t334 + t33 * t331) * t352 + ((t18 * t213 - t33 * t319 + (t213 * t78 - t87 * t319) * t358) * t310 + (t17 * t209 - t32 * t320 + (t209 * t77 - t86 * t320) * t360) * t311 + (t16 * t205 - t31 * t321 + (t205 * t76 - t85 * t321) * t362) * t312 + m(4) + (t205 ^ 2 + t209 ^ 2 + t213 ^ 2) * t141 * t255 ^ 2) * t229;];
tauX  = t13;
