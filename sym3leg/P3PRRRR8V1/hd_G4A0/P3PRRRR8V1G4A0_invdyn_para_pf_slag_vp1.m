% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:24:37
% EndTime: 2020-08-06 17:24:48
% DurationCPUTime: 10.64s
% Computational Cost: add. (44124->500), mult. (115560->1013), div. (7479->9), fcn. (142836->40), ass. (0->380)
t246 = cos(qJ(3,2));
t219 = 0.1e1 / t246;
t247 = cos(qJ(2,2));
t241 = sin(qJ(2,2));
t321 = t241 * t246;
t174 = pkin(2) * t321 - pkin(5) * t247;
t223 = sin(pkin(3));
t225 = cos(pkin(3));
t240 = sin(qJ(3,2));
t153 = t225 * t240 * pkin(2) + t174 * t223;
t421 = 0.1e1 / t153;
t357 = t421 * t219;
t242 = sin(qJ(3,1));
t248 = cos(qJ(3,1));
t432 = -rSges(3,1) * t248 + rSges(3,2) * t242;
t431 = -rSges(3,1) * t246 + rSges(3,2) * t240;
t238 = sin(qJ(3,3));
t244 = cos(qJ(3,3));
t430 = -rSges(3,1) * t244 + rSges(3,2) * t238;
t231 = legFrame(1,1);
t202 = sin(t231);
t208 = cos(t231);
t234 = legFrame(1,2);
t211 = sin(t234);
t214 = cos(t234);
t141 = g(1) * t211 + (-g(2) * t202 + g(3) * t208) * t214;
t228 = legFrame(1,3);
t199 = sin(t228);
t205 = cos(t228);
t344 = t208 * t211;
t347 = t202 * t211;
t398 = g(1) * t214;
t107 = -t199 * t398 + (-t199 * t347 + t205 * t208) * g(2) + (t199 * t344 + t202 * t205) * g(3);
t108 = t205 * t398 + (t199 * t208 + t205 * t347) * g(2) + (t199 * t202 - t205 * t344) * g(3);
t222 = sin(pkin(6));
t224 = cos(pkin(6));
t278 = t107 * t224 - t108 * t222;
t429 = -t141 * t225 + t278 * t223;
t230 = legFrame(2,1);
t201 = sin(t230);
t207 = cos(t230);
t233 = legFrame(2,2);
t210 = sin(t233);
t213 = cos(t233);
t140 = g(1) * t210 + (-g(2) * t201 + g(3) * t207) * t213;
t227 = legFrame(2,3);
t198 = sin(t227);
t204 = cos(t227);
t345 = t207 * t210;
t348 = t201 * t210;
t399 = g(1) * t213;
t105 = -t198 * t399 + (-t198 * t348 + t204 * t207) * g(2) + (t198 * t345 + t201 * t204) * g(3);
t106 = t204 * t399 + (t198 * t207 + t204 * t348) * g(2) + (t198 * t201 - t204 * t345) * g(3);
t280 = t105 * t224 - t106 * t222;
t428 = -t140 * t225 + t280 * t223;
t229 = legFrame(3,1);
t200 = sin(t229);
t206 = cos(t229);
t232 = legFrame(3,2);
t209 = sin(t232);
t212 = cos(t232);
t139 = g(1) * t209 + (-g(2) * t200 + g(3) * t206) * t212;
t226 = legFrame(3,3);
t197 = sin(t226);
t203 = cos(t226);
t346 = t206 * t209;
t349 = t200 * t209;
t400 = g(1) * t212;
t103 = -t197 * t400 + (-t197 * t349 + t203 * t206) * g(2) + (t197 * t346 + t200 * t203) * g(3);
t104 = t203 * t400 + (t197 * t206 + t203 * t349) * g(2) + (t197 * t200 - t203 * t346) * g(3);
t282 = t103 * t224 - t104 * t222;
t427 = -t139 * t225 + t282 * t223;
t426 = t141 * t223 + t278 * t225;
t425 = t140 * t223 + t280 * t225;
t424 = t139 * t223 + t282 * t225;
t245 = cos(qJ(2,3));
t239 = sin(qJ(2,3));
t330 = t225 * t239;
t164 = t222 * t245 + t224 * t330;
t167 = -t222 * t330 + t224 * t245;
t130 = t164 * t203 + t167 * t197;
t161 = -t197 * t222 + t203 * t224;
t337 = t223 * t244;
t100 = t130 * t238 + t161 * t337;
t322 = t239 * t244;
t173 = pkin(2) * t322 - pkin(5) * t245;
t331 = t225 * t238;
t152 = pkin(2) * t331 + t173 * t223;
t149 = 0.1e1 / t152;
t217 = 0.1e1 / t244;
t253 = xDP(3);
t254 = xDP(2);
t255 = xDP(1);
t343 = t212 * t255;
t158 = t197 * t224 + t203 * t222;
t110 = -t158 * t200 + t161 * t346;
t276 = t164 * t197 - t167 * t203;
t79 = (t130 * t346 - t276 * t200) * t238 + t110 * t337;
t115 = t158 * t206 + t161 * t349;
t82 = (-t130 * t349 - t276 * t206) * t238 - t115 * t337;
t55 = (-t100 * t343 + t253 * t79 + t254 * t82) * t217 * t149;
t52 = t55 ^ 2;
t329 = t225 * t241;
t165 = t222 * t247 + t224 * t329;
t168 = -t222 * t329 + t224 * t247;
t131 = t165 * t204 + t168 * t198;
t162 = -t198 * t222 + t204 * t224;
t335 = t223 * t246;
t101 = t131 * t240 + t162 * t335;
t342 = t213 * t255;
t159 = t198 * t224 + t204 * t222;
t112 = -t159 * t201 + t162 * t345;
t275 = t165 * t198 - t168 * t204;
t80 = (t131 * t345 - t275 * t201) * t240 + t112 * t335;
t116 = t159 * t207 + t162 * t348;
t83 = (-t131 * t348 - t275 * t207) * t240 - t116 * t335;
t56 = (-t101 * t342 + t253 * t80 + t254 * t83) * t357;
t53 = t56 ^ 2;
t249 = cos(qJ(2,1));
t243 = sin(qJ(2,1));
t327 = t225 * t243;
t166 = t222 * t249 + t224 * t327;
t169 = -t222 * t327 + t224 * t249;
t132 = t166 * t205 + t169 * t199;
t163 = -t199 * t222 + t205 * t224;
t333 = t223 * t248;
t102 = t132 * t242 + t163 * t333;
t320 = t243 * t248;
t175 = pkin(2) * t320 - pkin(5) * t249;
t328 = t225 * t242;
t154 = pkin(2) * t328 + t175 * t223;
t151 = 0.1e1 / t154;
t221 = 0.1e1 / t248;
t341 = t214 * t255;
t160 = t199 * t224 + t205 * t222;
t114 = -t160 * t202 + t163 * t344;
t274 = t166 * t199 - t169 * t205;
t81 = (t132 * t344 - t274 * t202) * t242 + t114 * t333;
t117 = t160 * t208 + t163 * t347;
t84 = (-t132 * t347 - t274 * t208) * t242 - t117 * t333;
t57 = (-t102 * t341 + t253 * t81 + t254 * t84) * t221 * t151;
t54 = t57 ^ 2;
t303 = t223 * t322;
t336 = t223 * t245;
t423 = 0.1e1 / (-pkin(5) * t336 + (t303 + t331) * pkin(2));
t301 = t223 * t320;
t332 = t223 * t249;
t422 = 0.1e1 / (-pkin(5) * t332 + (t301 + t328) * pkin(2));
t196 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t420 = 0.2e1 * t196;
t252 = m(2) * rSges(2,1);
t419 = m(3) * rSges(3,3);
t263 = 0.1e1 / pkin(2);
t358 = t423 * t217;
t306 = t263 * t358;
t109 = t158 * t346 + t161 * t200;
t326 = t225 * t245;
t403 = pkin(2) * t244;
t73 = (-t109 * t239 + t110 * t326) * t403 + pkin(5) * (t109 * t245 + t110 * t330);
t118 = t158 * t349 - t161 * t206;
t74 = -(t115 * t326 - t118 * t239) * t403 - pkin(5) * (t115 * t330 + t118 * t245);
t97 = (-t158 * t239 + t161 * t326) * t403 + pkin(5) * (t158 * t245 + t161 * t330);
t49 = (t253 * t73 + t254 * t74 - t97 * t343) * t306;
t418 = pkin(2) * t49;
t305 = t263 * t357;
t111 = t159 * t345 + t162 * t201;
t325 = t225 * t247;
t402 = pkin(2) * t246;
t75 = (-t111 * t241 + t112 * t325) * t402 + pkin(5) * (t111 * t247 + t112 * t329);
t119 = t159 * t348 - t162 * t207;
t76 = -(t116 * t325 - t119 * t241) * t402 - pkin(5) * (t116 * t329 + t119 * t247);
t98 = (-t159 * t241 + t162 * t325) * t402 + pkin(5) * (t159 * t247 + t162 * t329);
t50 = (t253 * t75 + t254 * t76 - t98 * t342) * t305;
t417 = pkin(2) * t50;
t356 = t422 * t221;
t304 = t263 * t356;
t113 = t160 * t344 + t163 * t202;
t324 = t225 * t249;
t401 = pkin(2) * t248;
t77 = (-t113 * t243 + t114 * t324) * t401 + pkin(5) * (t113 * t249 + t114 * t327);
t120 = t160 * t347 - t163 * t208;
t78 = -(t117 * t324 - t120 * t243) * t401 - pkin(5) * (t117 * t327 + t120 * t249);
t99 = (-t160 * t243 + t163 * t324) * t401 + pkin(5) * (t160 * t249 + t163 * t327);
t51 = (t253 * t77 + t254 * t78 - t99 * t341) * t304;
t416 = pkin(2) * t51;
t415 = pkin(5) * t55;
t414 = pkin(5) * t56;
t413 = pkin(5) * t57;
t259 = rSges(3,2) ^ 2;
t260 = rSges(3,1) ^ 2;
t186 = (-t259 + t260) * m(3) - Icges(3,1) + Icges(3,2);
t412 = t186 / 0.2e1;
t194 = rSges(3,2) * t419 - Icges(3,6);
t411 = -t194 / 0.4e1;
t195 = rSges(3,1) * t419 - Icges(3,5);
t410 = t195 / 0.4e1;
t390 = rSges(3,2) * t244;
t179 = rSges(3,1) * t238 + t390;
t133 = -t223 * t239 * t179 - t225 * t430;
t409 = m(3) * t133;
t389 = rSges(3,2) * t246;
t180 = rSges(3,1) * t240 + t389;
t134 = -t223 * t241 * t180 - t225 * t431;
t408 = m(3) * t134;
t388 = rSges(3,2) * t248;
t181 = rSges(3,1) * t242 + t388;
t135 = -t223 * t243 * t181 - t225 * t432;
t407 = m(3) * t135;
t216 = t244 ^ 2;
t406 = pkin(2) * t216;
t218 = t246 ^ 2;
t405 = pkin(2) * t218;
t220 = t248 ^ 2;
t404 = pkin(2) * t220;
t261 = pkin(5) ^ 2;
t262 = pkin(2) ^ 2;
t375 = t240 * t50;
t317 = pkin(2) * t375;
t397 = (-pkin(5) * t317 + (t218 * t262 + t261) * t56) * t56;
t387 = t423 * t55;
t386 = t422 * t57;
t385 = t212 * t97;
t384 = t213 * t98;
t383 = t214 * t99;
t382 = t217 * t79;
t381 = t217 * t82;
t380 = t219 * t80;
t379 = t219 * t83;
t378 = t221 * t81;
t377 = t221 * t84;
t376 = t238 * t49;
t374 = t242 * t51;
t373 = t100 * t212;
t372 = t101 * t213;
t371 = t102 * t214;
t256 = 0.2e1 * qJ(3,3);
t319 = t259 + t260;
t267 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t319) * m(3) / 0.2e1;
t121 = cos(t256) * t412 - t196 * sin(t256) + t267;
t370 = t121 * t217;
t257 = 0.2e1 * qJ(3,2);
t122 = cos(t257) * t412 - t196 * sin(t257) + t267;
t369 = t122 * t219;
t258 = 0.2e1 * qJ(3,1);
t123 = cos(t258) * t412 - t196 * sin(t258) + t267;
t368 = t123 * t221;
t170 = -m(3) * t430 + t252;
t193 = m(2) * rSges(2,2) - t419;
t136 = t170 * t245 - t193 * t239;
t367 = t136 * t223;
t171 = -m(3) * t431 + t252;
t137 = t171 * t247 - t193 * t241;
t366 = t137 * t223;
t172 = -m(3) * t432 + t252;
t138 = t172 * t249 - t193 * t243;
t365 = t138 * t223;
t155 = -t194 * t244 - t195 * t238;
t355 = t155 * t217;
t156 = -t194 * t246 - t195 * t240;
t354 = t156 * t219;
t157 = -t194 * t248 - t195 * t242;
t353 = t157 * t221;
t352 = t186 * t238;
t351 = t186 * t240;
t350 = t186 * t242;
t340 = t223 * t238;
t339 = t223 * t240;
t338 = t223 * t242;
t334 = t223 * t247;
t323 = t225 * t263;
t318 = pkin(2) * t376;
t316 = pkin(2) * t374;
t315 = t238 * t415;
t314 = t240 * t414;
t313 = t242 * t413;
t312 = t217 * t373;
t311 = t219 * t372;
t310 = t221 * t371;
t309 = t217 * t367;
t308 = t219 * t366;
t307 = t221 * t365;
t302 = t223 * t321;
t300 = t155 * t306;
t191 = t319 * m(3) + Icges(3,3);
t299 = t191 * t306;
t298 = t156 * t305;
t297 = t191 * t305;
t296 = t157 * t304;
t295 = t191 * t304;
t294 = t306 * t409;
t293 = t305 * t408;
t292 = t304 * t407;
t291 = t306 * t385;
t290 = t305 * t384;
t289 = t304 * t383;
t281 = t103 * t222 + t104 * t224;
t279 = t105 * t222 + t106 * t224;
t277 = t107 * t222 + t108 * t224;
t273 = pkin(2) * t340 - t173 * t225;
t272 = pkin(2) * t339 - t174 * t225;
t271 = pkin(2) * t338 - t175 * t225;
t235 = xDDP(3);
t236 = xDDP(2);
t237 = xDDP(1);
t270 = t235 * t73 + t236 * t74 - t237 * t385;
t269 = t235 * t75 + t236 * t76 - t237 * t384;
t268 = t235 * t77 + t236 * t78 - t237 * t383;
t266 = t424 * t239 + t281 * t245;
t265 = t425 * t241 + t279 * t247;
t264 = t426 * t243 + t277 * t249;
t215 = m(1) + m(2) + m(3);
t178 = pkin(5) * t243 + t249 * t401;
t177 = pkin(5) * t241 + t247 * t402;
t176 = pkin(5) * t239 + t245 * t403;
t129 = -t178 * t222 + t271 * t224;
t128 = -t177 * t222 + t272 * t224;
t127 = -t176 * t222 + t273 * t224;
t126 = t178 * t224 + t271 * t222;
t125 = t177 * t224 + t272 * t222;
t124 = t176 * t224 + t273 * t222;
t96 = -t126 * t199 + t129 * t205;
t95 = -t125 * t198 + t128 * t204;
t94 = -t124 * t197 + t127 * t203;
t93 = (t271 * t160 + t163 * t178) * t214 + t211 * t154;
t92 = (t272 * t159 + t162 * t177) * t213 + t210 * t153;
t91 = (t273 * t158 + t161 * t176) * t212 + t209 * t152;
t87 = -t154 * t214 + (t126 * t205 + t129 * t199) * t211;
t86 = -t153 * t213 + (t125 * t204 + t128 * t198) * t210;
t85 = -t152 * t212 + (t124 * t203 + t127 * t197) * t209;
t72 = -t202 * t96 - t208 * t87;
t71 = t202 * t87 - t208 * t96;
t70 = -t201 * t95 - t207 * t86;
t69 = t201 * t86 - t207 * t95;
t68 = -t200 * t94 - t206 * t85;
t67 = t200 * t85 - t206 * t94;
t63 = -t289 * t407 + (t215 * t93 - t307 * t371) * t151;
t62 = -t290 * t408 + (t215 * t92 - t308 * t372) * t421;
t61 = -t291 * t409 + (t215 * t91 - t309 * t373) * t149;
t60 = -t157 * t289 + (-t123 * t310 + t93 * t365) * t151;
t59 = -t156 * t290 + (-t122 * t311 + t92 * t366) * t421;
t58 = -t155 * t291 + (-t121 * t312 + t91 * t367) * t149;
t48 = t51 ^ 2;
t47 = t50 ^ 2;
t46 = t49 ^ 2;
t39 = t77 * t292 + (t215 * t72 + t81 * t307) * t151;
t38 = t78 * t292 + (t215 * t71 + t84 * t307) * t151;
t37 = t75 * t293 + (t215 * t70 + t80 * t308) * t421;
t36 = t76 * t293 + (t215 * t69 + t83 * t308) * t421;
t35 = t73 * t294 + (t215 * t68 + t79 * t309) * t149;
t34 = t74 * t294 + (t215 * t67 + t82 * t309) * t149;
t33 = t77 * t296 + (t72 * t365 + t81 * t368) * t151;
t32 = t78 * t296 + (t71 * t365 + t84 * t368) * t151;
t31 = t75 * t298 + (t70 * t366 + t80 * t369) * t421;
t30 = t76 * t298 + (t69 * t366 + t83 * t369) * t421;
t29 = t73 * t300 + (t68 * t367 + t79 * t370) * t149;
t28 = t74 * t300 + (t67 * t367 + t82 * t370) * t149;
t24 = t313 - t416;
t23 = t314 - t417;
t22 = t315 - t418;
t21 = -pkin(5) * t316 + (t220 * t262 + t261) * t57;
t19 = -pkin(5) * t318 + (t216 * t262 + t261) * t55;
t18 = (t21 * t57 - t24 * t416) * t422;
t17 = (t23 * t417 - t397) * t421;
t16 = (t19 * t55 - t22 * t418) * t423;
t15 = (t21 * t323 * t386 + (-t51 * t175 * t338 + t225 * (t51 * t404 - t313)) * t151 * t51) * t221;
t14 = (t323 * t397 + (-t50 * t174 * t339 + t225 * (t50 * t405 - t314)) * t50) * t357;
t13 = (t19 * t323 * t387 + (-t49 * t173 * t340 + t225 * (t49 * t406 - t315)) * t149 * t49) * t217;
t12 = (((t225 * t51 + t57 * t332) * t404 - (t316 - t413) * t301 + t225 * t24) * t386 + (t51 * t332 + (t220 * t225 - t242 * t301 - t225) * t57) * t422 * t416) * t221;
t11 = (((t225 * t50 + t56 * t334) * t405 - (t317 - t414) * t302 + t225 * t23) * t56 + (t50 * t334 + (t218 * t225 - t240 * t302 - t225) * t56) * t417) * t357;
t10 = (((t225 * t49 + t55 * t336) * t406 - (t318 - t415) * t303 + t225 * t22) * t387 + (t49 * t336 + (t216 * t225 - t238 * t303 - t225) * t55) * t423 * t418) * t217;
t9 = t18 * t407 - t157 * t12 - t191 * t15 + t54 * (t220 * t420 + t248 * t350 - t196) + m(3) * ((t429 * rSges(3,1) + t264 * rSges(3,2)) * t248 + t242 * (t264 * rSges(3,1) - t429 * rSges(3,2)));
t8 = -t17 * t408 - t156 * t11 - t191 * t14 + t53 * (t218 * t420 + t246 * t351 - t196) + m(3) * ((t428 * rSges(3,1) + t265 * rSges(3,2)) * t246 + t240 * (t265 * rSges(3,1) - t428 * rSges(3,2)));
t7 = t16 * t409 - t155 * t10 - t191 * t13 + t52 * (t216 * t420 + t244 * t352 - t196) + m(3) * ((t427 * rSges(3,1) + t266 * rSges(3,2)) * t244 + t238 * (t266 * rSges(3,1) - t427 * rSges(3,2)));
t6 = t18 * t365 - t123 * t12 - t157 * t15 - 0.4e1 * ((t57 * t350 / 0.2e1 + t51 * t410) * t248 + t374 * t411 + (t220 - 0.1e1 / 0.2e1) * t57 * t196) * t51 + (-t172 * t426 + t277 * t193) * t249 - t243 * (-t277 * t172 - t193 * t426);
t5 = -t17 * t366 - t122 * t11 - t156 * t14 - 0.4e1 * ((t56 * t351 / 0.2e1 + t50 * t410) * t246 + t375 * t411 + (t218 - 0.1e1 / 0.2e1) * t56 * t196) * t50 + (-t171 * t425 + t279 * t193) * t247 - t241 * (-t279 * t171 - t193 * t425);
t4 = t16 * t367 - t121 * t10 - t155 * t13 - 0.4e1 * ((t55 * t352 / 0.2e1 + t49 * t410) * t244 + t376 * t411 + (t216 - 0.1e1 / 0.2e1) * t55 * t196) * t49 + (-t170 * t424 + t281 * t193) * t245 - t239 * (-t281 * t170 - t193 * t424);
t3 = (-t138 * t12 + (-t193 * t249 - t243 * t252) * t54) * t223 + (t18 - t141) * t215 + (-t135 * t15 + (-0.2e1 * t57 * t249 * (rSges(3,1) * t374 + t51 * t388) + t432 * t243 * (t54 + t48)) * t223 - t48 * t225 * t181) * m(3);
t2 = (-t137 * t11 + (-t193 * t247 - t241 * t252) * t53) * t223 + (-t17 - t140) * t215 + (-t134 * t14 + (-0.2e1 * t56 * t247 * (rSges(3,1) * t375 + t50 * t389) + t431 * t241 * (t53 + t47)) * t223 - t47 * t225 * t180) * m(3);
t1 = (-t136 * t10 + (-t193 * t245 - t239 * t252) * t52) * t223 + (t16 - t139) * t215 + (-t133 * t13 + (-0.2e1 * t55 * t245 * (rSges(3,1) * t376 + t49 * t390) + t430 * t239 * (t52 + t46)) * t223 - t46 * t225 * t179) * m(3);
t20 = [(-g(1) + t237) * m(4) + ((-t60 * t310 + t63 * t93) * t237 + (t60 * t377 + t63 * t71) * t236 + (t60 * t378 + t63 * t72) * t235 + t93 * t3 - t6 * t310) * t151 + ((-t59 * t311 + t62 * t92) * t237 + (t59 * t379 + t62 * t69) * t236 + (t59 * t380 + t62 * t70) * t235 + t92 * t2 - t5 * t311) * t421 + ((-t9 * t383 + t268 * (-t191 * t289 + (-t157 * t310 + t93 * t407) * t151)) * t356 + (-t8 * t384 + t269 * (-t191 * t290 + (-t156 * t311 + t92 * t408) * t421)) * t357 + (-t7 * t385 + t270 * (-t191 * t291 + (-t155 * t312 + t91 * t409) * t149)) * t358) * t263 + ((-t58 * t312 + t61 * t91) * t237 + (t58 * t381 + t61 * t67) * t236 + (t58 * t382 + t61 * t68) * t235 + t91 * t1 - t4 * t312) * t149; (-g(2) + t236) * m(4) + ((-t32 * t310 + t38 * t93) * t237 + (t32 * t377 + t38 * t71) * t236 + (t32 * t378 + t38 * t72) * t235 + t71 * t3 + t6 * t377) * t151 + ((-t30 * t311 + t36 * t92) * t237 + (t30 * t379 + t36 * t69) * t236 + (t30 * t380 + t36 * t70) * t235 + t69 * t2 + t5 * t379) * t421 + ((t78 * t9 + t268 * (t78 * t295 + (t84 * t353 + t71 * t407) * t151)) * t356 + (t76 * t8 + t269 * (t76 * t297 + (t83 * t354 + t69 * t408) * t421)) * t357 + (t7 * t74 + t270 * (t74 * t299 + (t82 * t355 + t67 * t409) * t149)) * t358) * t263 + ((-t28 * t312 + t34 * t91) * t237 + (t28 * t381 + t34 * t67) * t236 + (t28 * t382 + t34 * t68) * t235 + t67 * t1 + t4 * t381) * t149; (-g(3) + t235) * m(4) + ((-t33 * t310 + t39 * t93) * t237 + (t33 * t377 + t39 * t71) * t236 + (t33 * t378 + t39 * t72) * t235 + t72 * t3 + t6 * t378) * t151 + ((-t31 * t311 + t37 * t92) * t237 + (t31 * t379 + t37 * t69) * t236 + (t31 * t380 + t37 * t70) * t235 + t70 * t2 + t5 * t380) * t421 + ((t77 * t9 + t268 * (t77 * t295 + (t81 * t353 + t72 * t407) * t151)) * t356 + (t75 * t8 + t269 * (t75 * t297 + (t80 * t354 + t70 * t408) * t421)) * t357 + (t7 * t73 + t270 * (t73 * t299 + (t79 * t355 + t68 * t409) * t149)) * t358) * t263 + ((-t29 * t312 + t35 * t91) * t237 + (t29 * t381 + t35 * t67) * t236 + (t29 * t382 + t35 * t68) * t235 + t68 * t1 + t4 * t382) * t149;];
tauX  = t20;
