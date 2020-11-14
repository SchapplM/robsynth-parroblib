% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
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
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:22
% EndTime: 2020-08-07 10:58:28
% DurationCPUTime: 5.80s
% Computational Cost: add. (7649->497), mult. (16307->946), div. (5940->13), fcn. (13474->34), ass. (0->374)
t236 = xP(4);
t183 = sin(t236);
t184 = cos(t236);
t246 = koppelP(1,2);
t250 = koppelP(1,1);
t142 = t183 * t250 + t184 * t246;
t146 = -t183 * t246 + t184 * t250;
t232 = xDP(4);
t202 = t232 ^ 2;
t211 = xDDP(4);
t214 = xDDP(1);
t112 = -t142 * t211 - t146 * t202 + t214;
t210 = legFrame(1,2);
t182 = cos(t210);
t376 = t112 * t182;
t213 = xDDP(2);
t108 = -t142 * t202 + t146 * t211 + t213;
t178 = sin(t210);
t384 = t108 * t178;
t417 = t376 - t384;
t245 = koppelP(2,2);
t249 = koppelP(2,1);
t141 = t183 * t249 + t184 * t245;
t145 = -t183 * t245 + t184 * t249;
t111 = -t141 * t211 - t145 * t202 + t214;
t209 = legFrame(2,2);
t181 = cos(t209);
t378 = t111 * t181;
t107 = -t141 * t202 + t145 * t211 + t213;
t177 = sin(t209);
t386 = t107 * t177;
t416 = t378 - t386;
t244 = koppelP(3,2);
t248 = koppelP(3,1);
t140 = t183 * t248 + t184 * t244;
t144 = -t183 * t244 + t184 * t248;
t110 = -t140 * t211 - t144 * t202 + t214;
t208 = legFrame(3,2);
t180 = cos(t208);
t380 = t110 * t180;
t106 = -t140 * t202 + t144 * t211 + t213;
t176 = sin(t208);
t388 = t106 * t176;
t415 = t380 - t388;
t243 = koppelP(4,2);
t247 = koppelP(4,1);
t139 = t183 * t247 + t184 * t243;
t143 = -t183 * t243 + t184 * t247;
t109 = -t139 * t211 - t143 * t202 + t214;
t207 = legFrame(4,2);
t179 = cos(t207);
t382 = t109 * t179;
t105 = -t139 * t202 + t143 * t211 + t213;
t175 = sin(t207);
t390 = t105 * t175;
t414 = t382 - t390;
t204 = sin(qJ(2,4));
t185 = 0.1e1 / t204;
t205 = cos(qJ(3,4));
t187 = 0.1e1 / t205;
t233 = xDP(3);
t253 = 0.1e1 / pkin(2);
t234 = xDP(2);
t235 = xDP(1);
t298 = (t143 * t232 + t234) * t175 - (-t139 * t232 + t235) * t179;
t188 = 0.1e1 / t205 ^ 2;
t203 = sin(qJ(3,4));
t206 = cos(qJ(2,4));
t317 = t188 * t203 * t206;
t60 = (-t187 * t233 - t298 * t317) * t253 * t185;
t59 = t60 ^ 2;
t216 = sin(qJ(2,3));
t190 = 0.1e1 / t216;
t221 = cos(qJ(3,3));
t194 = 0.1e1 / t221;
t297 = (t144 * t232 + t234) * t176 - (-t140 * t232 + t235) * t180;
t195 = 0.1e1 / t221 ^ 2;
t215 = sin(qJ(3,3));
t222 = cos(qJ(2,3));
t313 = t195 * t215 * t222;
t70 = (-t194 * t233 - t297 * t313) * t253 * t190;
t61 = t70 ^ 2;
t218 = sin(qJ(2,2));
t191 = 0.1e1 / t218;
t223 = cos(qJ(3,2));
t197 = 0.1e1 / t223;
t296 = (t145 * t232 + t234) * t177 - (-t141 * t232 + t235) * t181;
t198 = 0.1e1 / t223 ^ 2;
t217 = sin(qJ(3,2));
t224 = cos(qJ(2,2));
t312 = t198 * t217 * t224;
t71 = (-t197 * t233 - t296 * t312) * t253 * t191;
t62 = t71 ^ 2;
t220 = sin(qJ(2,1));
t192 = 0.1e1 / t220;
t225 = cos(qJ(3,1));
t200 = 0.1e1 / t225;
t295 = (t146 * t232 + t234) * t178 - (-t142 * t232 + t235) * t182;
t201 = 0.1e1 / t225 ^ 2;
t219 = sin(qJ(3,1));
t226 = cos(qJ(2,1));
t311 = t201 * t219 * t226;
t72 = (-t200 * t233 - t295 * t311) * t253 * t192;
t63 = t72 ^ 2;
t174 = (m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4));
t413 = 2 * t174;
t231 = m(2) * rSges(2,1);
t412 = rSges(3,3) * m(3);
t411 = g(3) * rSges(3,3);
t251 = rSges(3,2) ^ 2;
t252 = rSges(3,1) ^ 2;
t165 = (-t251 + t252) * m(3) - Icges(3,1) + Icges(3,2);
t410 = t165 / 0.2e1;
t172 = rSges(3,2) * t412 - Icges(3,6);
t409 = -t172 / 0.4e1;
t173 = rSges(3,1) * t412 - Icges(3,5);
t408 = t173 / 0.4e1;
t407 = rSges(3,1) * t205;
t406 = rSges(3,1) * t221;
t405 = rSges(3,1) * t223;
t404 = rSges(3,1) * t225;
t353 = t187 * t253;
t85 = t298 * t353;
t81 = t85 ^ 2;
t403 = t187 * t81;
t346 = t194 * t253;
t86 = t297 * t346;
t82 = t86 ^ 2;
t402 = t194 * t82;
t345 = t197 * t253;
t87 = t296 * t345;
t83 = t87 ^ 2;
t401 = t197 * t83;
t344 = t200 * t253;
t88 = t295 * t344;
t84 = t88 ^ 2;
t400 = t200 * t84;
t399 = t203 * t85;
t398 = t206 * t85;
t397 = t215 * t86;
t396 = t217 * t87;
t395 = t219 * t88;
t394 = t222 * t86;
t393 = t224 * t87;
t392 = t226 * t88;
t343 = t204 * t205;
t126 = t175 * t203 + t179 * t343;
t391 = t105 * t126;
t338 = t216 * t221;
t134 = t176 * t215 + t180 * t338;
t389 = t106 * t134;
t337 = t218 * t223;
t135 = t177 * t217 + t181 * t337;
t387 = t107 * t135;
t336 = t220 * t225;
t136 = t178 * t219 + t182 * t336;
t385 = t108 * t136;
t125 = t175 * t343 - t179 * t203;
t383 = t109 * t125;
t131 = t176 * t338 - t180 * t215;
t381 = t110 * t131;
t132 = t177 * t337 - t181 * t217;
t379 = t111 * t132;
t133 = t178 * t336 - t182 * t219;
t377 = t112 * t133;
t155 = rSges(3,1) * t203 + rSges(3,2) * t205;
t375 = t155 * t187;
t374 = t155 * t204;
t156 = rSges(3,1) * t215 + rSges(3,2) * t221;
t373 = t156 * t194;
t372 = t156 * t216;
t157 = rSges(3,1) * t217 + rSges(3,2) * t223;
t371 = t157 * t197;
t370 = t157 * t218;
t158 = rSges(3,1) * t219 + rSges(3,2) * t225;
t369 = t158 * t200;
t368 = t158 * t220;
t367 = t165 * t203;
t366 = t165 * t215;
t365 = t165 * t217;
t364 = t165 * t219;
t363 = t175 * t253;
t362 = t176 * t253;
t361 = t177 * t253;
t360 = t178 * t253;
t359 = t179 * t253;
t358 = t180 * t253;
t357 = t181 * t253;
t356 = t182 * t253;
t355 = t185 * t187;
t354 = t185 * t206;
t352 = t190 * t194;
t351 = t190 * t222;
t350 = t191 * t197;
t349 = t191 * t224;
t348 = t192 * t200;
t347 = t192 * t226;
t212 = xDDP(3);
t342 = t206 * t212;
t341 = t212 * t222;
t340 = t212 * t224;
t339 = t212 * t226;
t335 = (t251 + t252);
t334 = m(3) * t375;
t333 = m(3) * t374;
t332 = m(3) * t373;
t331 = m(3) * t372;
t330 = m(3) * t371;
t329 = m(3) * t370;
t328 = m(3) * t369;
t327 = m(3) * t368;
t326 = t203 * t403;
t325 = t215 * t402;
t324 = t217 * t401;
t323 = t219 * t400;
t171 = m(2) * rSges(2,2) - t412;
t306 = -rSges(3,2) * t203 + t407;
t113 = (t306 * m(3) + t231) * t206 - t204 * t171;
t322 = t113 * t355;
t305 = -rSges(3,2) * t215 + t406;
t114 = (t305 * m(3) + t231) * t222 - t216 * t171;
t321 = t114 * t352;
t304 = -rSges(3,2) * t217 + t405;
t115 = (t304 * m(3) + t231) * t224 - t218 * t171;
t320 = t115 * t350;
t303 = -rSges(3,2) * t219 + t404;
t116 = (t303 * m(3) + t231) * t226 - t220 * t171;
t319 = t116 * t348;
t189 = m(1) + m(2) + m(3);
t318 = t189 * t355;
t316 = t189 * t352;
t315 = t189 * t350;
t314 = t189 * t348;
t310 = t185 * t317;
t309 = t190 * t313;
t308 = t191 * t312;
t307 = t192 * t311;
t147 = g(1) * t175 + g(2) * t179;
t302 = g(3) * t206 + t147 * t204;
t148 = g(1) * t176 + g(2) * t180;
t301 = g(3) * t222 + t148 * t216;
t149 = g(1) * t177 + g(2) * t181;
t300 = g(3) * t224 + t149 * t218;
t150 = g(1) * t178 + g(2) * t182;
t299 = g(3) * t226 + t150 * t220;
t294 = t139 * t179 + t143 * t175;
t293 = t140 * t180 + t144 * t176;
t292 = t141 * t181 + t145 * t177;
t291 = t142 * t182 + t146 * t178;
t290 = t187 * (t383 + t391);
t289 = t194 * (t381 + t389);
t288 = t197 * (t379 + t387);
t287 = (t377 + t385) * t200;
t127 = -t172 * t205 - t173 * t203;
t13 = ((t205 * t206 * t60 - t204 * t399) * t187 * t60 + (-t203 * t60 * t343 + t398) * t188 * t85) * t185;
t186 = t205 ^ 2;
t229 = rSges(2,2) * g(3);
t230 = rSges(2,1) * g(3);
t29 = (-t205 * t59 - t403) * t185 * pkin(2);
t239 = 0.2e1 * qJ(3,4);
t254 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (2 * rSges(3,3) ^ 2 + t335) * m(3) / 0.2e1;
t89 = cos(t239) * t410 - t174 * sin(t239) + t254;
t5 = -t113 * t29 - t89 * t13 + t127 * t326 - 0.4e1 * ((t60 * t367 / 0.2e1 + t85 * t408) * t205 + t399 * t409 + (t186 - 0.1e1 / 0.2e1) * t60 * t174) * t85 + (m(2) * (-rSges(2,1) * t147 + t229) + (-t306 * t147 - t411) * m(3)) * t206 + t204 * (m(2) * (rSges(2,2) * t147 + t230) + (-t147 * rSges(3,3) + t306 * g(3)) * m(3));
t151 = g(1) * t179 - g(2) * t175;
t170 = t335 * m(3) + Icges(3,3);
t9 = -t127 * t13 + t170 * t326 + t59 * (t186 * t413 + t205 * t367 + Icges(3,4)) + (t29 * t374 + (-t203 * t151 + t302 * t205) * rSges(3,2) + (-t59 * rSges(3,2) + t151 * t205 + t203 * t302) * rSges(3,1)) * m(3);
t286 = -t187 * t9 + t5 * t310;
t128 = -t172 * t221 - t173 * t215;
t14 = ((t221 * t222 * t70 - t216 * t397) * t194 * t70 + (-t215 * t70 * t338 + t394) * t195 * t86) * t190;
t152 = g(1) * t180 - g(2) * t176;
t193 = t221 ^ 2;
t30 = (-t221 * t61 - t402) * t190 * pkin(2);
t10 = -t128 * t14 + t170 * t325 + t61 * (t193 * t413 + t221 * t366 + Icges(3,4)) + (t30 * t372 + (-t215 * t152 + t301 * t221) * rSges(3,2) + (-t61 * rSges(3,2) + t152 * t221 + t215 * t301) * rSges(3,1)) * m(3);
t240 = 0.2e1 * qJ(3,3);
t90 = cos(t240) * t410 - t174 * sin(t240) + t254;
t6 = -t114 * t30 - t90 * t14 + t128 * t325 - 0.4e1 * ((t70 * t366 / 0.2e1 + t86 * t408) * t221 + t397 * t409 + (t193 - 0.1e1 / 0.2e1) * t70 * t174) * t86 + (m(2) * (-rSges(2,1) * t148 + t229) + (-t305 * t148 - t411) * m(3)) * t222 + t216 * (m(2) * (rSges(2,2) * t148 + t230) + (-t148 * rSges(3,3) + t305 * g(3)) * m(3));
t285 = -t194 * t10 + t6 * t309;
t129 = -t172 * t223 - t173 * t217;
t15 = ((t223 * t224 * t71 - t218 * t396) * t197 * t71 + (-t217 * t71 * t337 + t393) * t198 * t87) * t191;
t153 = g(1) * t181 - g(2) * t177;
t196 = t223 ^ 2;
t31 = (-t223 * t62 - t401) * t191 * pkin(2);
t11 = -t129 * t15 + t170 * t324 + t62 * (t196 * t413 + t223 * t365 + Icges(3,4)) + (t31 * t370 + (-t217 * t153 + t300 * t223) * rSges(3,2) + (-t62 * rSges(3,2) + t153 * t223 + t217 * t300) * rSges(3,1)) * m(3);
t241 = 0.2e1 * qJ(3,2);
t91 = cos(t241) * t410 - t174 * sin(t241) + t254;
t7 = -t115 * t31 - t91 * t15 + t129 * t324 - 0.4e1 * ((t71 * t365 / 0.2e1 + t87 * t408) * t223 + t396 * t409 + (t196 - 0.1e1 / 0.2e1) * t71 * t174) * t87 + (m(2) * (-rSges(2,1) * t149 + t229) + (-t304 * t149 - t411) * m(3)) * t224 + t218 * (m(2) * (rSges(2,2) * t149 + t230) + (-t149 * rSges(3,3) + t304 * g(3)) * m(3));
t284 = -t197 * t11 + t7 * t308;
t130 = -t172 * t225 - t173 * t219;
t154 = g(1) * t182 - g(2) * t178;
t16 = ((t225 * t226 * t72 - t220 * t395) * t200 * t72 + (-t219 * t72 * t336 + t392) * t201 * t88) * t192;
t199 = t225 ^ 2;
t32 = (-t225 * t63 - t400) * t192 * pkin(2);
t12 = -t130 * t16 + t170 * t323 + t63 * (t199 * t413 + t225 * t364 + Icges(3,4)) + (t32 * t368 + (-t219 * t154 + t299 * t225) * rSges(3,2) + (-t63 * rSges(3,2) + t154 * t225 + t219 * t299) * rSges(3,1)) * m(3);
t242 = 0.2e1 * qJ(3,1);
t92 = cos(t242) * t410 - t174 * sin(t242) + t254;
t8 = -t116 * t32 - t92 * t16 + t130 * t323 - 0.4e1 * ((t72 * t364 / 0.2e1 + t88 * t408) * t225 + t395 * t409 + (t199 - 0.1e1 / 0.2e1) * t72 * t174) * t88 + (m(2) * (-rSges(2,1) * t150 + t229) + (-t303 * t150 - t411) * m(3)) * t226 + t220 * (m(2) * (rSges(2,2) * t150 + t230) + (-t150 * rSges(3,3) + t303 * g(3)) * m(3));
t283 = -t200 * t12 + t8 * t307;
t262 = t127 * t310 - t170 * t187;
t266 = -t127 * t187 + t89 * t310;
t33 = t125 * t322 + t266 * t359;
t280 = -t187 * (-t125 * t334 + t262 * t359) + t33 * t310;
t34 = t126 * t322 - t266 * t363;
t279 = -t187 * (-t126 * t334 - t262 * t363) + t34 * t310;
t261 = t128 * t309 - t170 * t194;
t265 = -t128 * t194 + t90 * t309;
t36 = t131 * t321 + t265 * t358;
t277 = -t194 * (-t131 * t332 + t261 * t358) + t36 * t309;
t39 = t134 * t321 - t265 * t362;
t276 = -t194 * (-t134 * t332 - t261 * t362) + t39 * t309;
t260 = t129 * t308 - t170 * t197;
t264 = -t129 * t197 + t91 * t308;
t37 = t132 * t320 + t264 * t357;
t275 = -t197 * (-t132 * t330 + t260 * t357) + t37 * t308;
t40 = t135 * t320 - t264 * t361;
t274 = -t197 * (-t135 * t330 - t260 * t361) + t40 * t308;
t259 = t130 * t307 - t170 * t200;
t263 = -t130 * t200 + t92 * t307;
t38 = t133 * t319 + t263 * t356;
t272 = -t200 * (-t133 * t328 + t259 * t356) + t38 * t307;
t41 = t136 * t319 - t263 * t360;
t271 = -t200 * (-t136 * t328 - t259 * t360) + t41 * t307;
t258 = t113 * t310 + t187 * t333;
t257 = t114 * t309 + t194 * t331;
t256 = t115 * t308 + t197 * t329;
t255 = t116 * t307 + t200 * t327;
t238 = rSges(4,1);
t237 = rSges(4,2);
t138 = -t183 * t237 + t184 * t238;
t137 = t183 * t238 + t184 * t237;
t100 = (-t116 * t344 + t189 * t226) * t192;
t99 = (-t115 * t345 + t189 * t224) * t191;
t98 = (-t114 * t346 + t189 * t222) * t190;
t97 = t291 * t344;
t96 = t292 * t345;
t95 = t293 * t346;
t94 = t294 * t353;
t93 = (-t113 * t353 + t189 * t206) * t185;
t80 = t291 * t253 * t307;
t79 = t292 * t253 * t308;
t78 = t293 * t253 * t309;
t77 = t294 * t253 * t310;
t76 = (-t133 * t142 + t136 * t146) * t348;
t75 = (-t132 * t141 + t135 * t145) * t350;
t74 = (-t131 * t140 + t134 * t144) * t352;
t73 = (-t125 * t139 + t126 * t143) * t355;
t56 = (t116 * t226 - t92 * t344) * t192;
t55 = (t115 * t224 - t91 * t345) * t191;
t54 = (t114 * t222 - t90 * t346) * t190;
t53 = (t113 * t206 - t89 * t353) * t185;
t52 = t136 * t314 - t255 * t360;
t51 = t135 * t315 - t256 * t361;
t50 = t134 * t316 - t257 * t362;
t49 = t133 * t314 + t255 * t356;
t48 = t132 * t315 + t256 * t357;
t47 = t131 * t316 + t257 * t358;
t46 = t126 * t318 - t258 * t363;
t45 = t125 * t318 + t258 * t359;
t44 = t63 + t84;
t43 = t62 + t83;
t42 = t61 + t82;
t35 = t59 + t81;
t25 = -t116 * t80 + t189 * t76 - t97 * t327;
t24 = -t115 * t79 + t189 * t75 - t96 * t329;
t23 = -t114 * t78 + t189 * t74 - t95 * t331;
t21 = -t113 * t77 + t189 * t73 - t94 * t333;
t20 = t116 * t76 + t130 * t97 - t80 * t92;
t19 = t115 * t75 + t129 * t96 - t79 * t91;
t18 = t114 * t74 + t128 * t95 - t78 * t90;
t17 = t113 * t73 + t127 * t94 - t77 * t89;
t4 = -t116 * t16 + (-t32 - t150) * t189 + (-0.2e1 * t72 * t158 * t392 + (-t44 * t404 + (rSges(3,2) * t44 - t84 * t369) * t219) * t220) * m(3) + (-t226 * t171 - t220 * t231) * t63;
t3 = -t115 * t15 + (-t31 - t149) * t189 + (-0.2e1 * t71 * t157 * t393 + (-t43 * t405 + (rSges(3,2) * t43 - t83 * t371) * t217) * t218) * m(3) + (-t224 * t171 - t218 * t231) * t62;
t2 = -t114 * t14 + (-t30 - t148) * t189 + (-0.2e1 * t70 * t156 * t394 + (-t42 * t406 + (rSges(3,2) * t42 - t82 * t373) * t215) * t216) * m(3) + (-t222 * t171 - t216 * t231) * t61;
t1 = -t113 * t13 + (-t29 - t147) * t189 + (-0.2e1 * t60 * t155 * t398 + (-t35 * t407 + (rSges(3,2) * t35 - t81 * t375) * t203) * t204) * m(3) + (-t206 * t171 - t204 * t231) * t59;
t22 = [(t49 * t339 + (t49 * t385 + (t112 * t49 + t4) * t133) * t200) * t192 + (t48 * t340 + (t48 * t387 + (t111 * t48 + t3) * t132) * t197) * t191 + (t47 * t341 + (t47 * t389 + (t110 * t47 + t2) * t131) * t194) * t190 + (t45 * t342 + (t45 * t391 + (t109 * t45 + t1) * t125) * t187) * t185 + (-t137 * t211 - t202 * t138 - g(1) + t214) * m(4) + (-t272 * t384 - t275 * t386 - t277 * t388 - t280 * t390 + (-t33 * t355 - t38 * t348 - t37 * t350 - t36 * t352) * t212 + (t272 * t112 + t283) * t182 + (t275 * t111 + t284) * t181 + (t277 * t110 + t285) * t180 + (t280 * t109 + t286) * t179) * t253; (t52 * t339 + (t52 * t377 + (t108 * t52 + t4) * t136) * t200) * t192 + (t51 * t340 + (t51 * t379 + (t107 * t51 + t3) * t135) * t197) * t191 + (t50 * t341 + (t50 * t381 + (t106 * t50 + t2) * t134) * t194) * t190 + (t46 * t342 + (t46 * t383 + (t105 * t46 + t1) * t126) * t187) * t185 + (-t137 * t202 + t138 * t211 - g(2) + t213) * m(4) + (t271 * t376 + t274 * t378 + t276 * t380 + t279 * t382 + (-t34 * t355 - t348 * t41 - t350 * t40 - t352 * t39) * t212 + (-t108 * t271 - t283) * t178 + (-t107 * t274 - t284) * t177 + (-t106 * t276 - t285) * t176 + (-t105 * t279 - t286) * t175) * t253; -m(4) * g(3) + (t100 * t287 + t226 * t4) * t192 + (t224 * t3 + t288 * t99) * t191 + (t222 * t2 + t289 * t98) * t190 + (t206 * t1 + t290 * t93) * t185 + (t100 * t347 + t349 * t99 + t351 * t98 + t93 * t354 + m(4)) * t212 + (-t5 * t355 - t6 * t352 - t7 * t350 - t8 * t348 + (-t348 * t56 - t350 * t55 - t352 * t54 - t53 * t355) * t212 + t414 * (-(-t206 * m(3) * t155 - t127 * t185 * t353) * t187 + t53 * t310) + t415 * (-(-t222 * m(3) * t156 - t128 * t190 * t346) * t194 + t54 * t309) + t416 * (-(-t224 * m(3) * t157 - t129 * t191 * t345) * t197 + t55 * t308) + t417 * (-(-t226 * m(3) * t158 - t130 * t192 * t344) * t200 + t56 * t307)) * t253; Icges(4,3) * t211 + t73 * t1 + t95 * t10 + t96 * t11 + t97 * t12 + t74 * t2 + t75 * t3 + t76 * t4 - t77 * t5 - t78 * t6 - t79 * t7 - t80 * t8 + t94 * t9 + t25 * t192 * t287 + t24 * t191 * t288 + t23 * t190 * t289 + t21 * t185 * t290 + ((t237 ^ 2 + t238 ^ 2) * t211 + (g(1) * t237 - g(2) * t238) * t184 + t183 * (g(1) * t238 + g(2) * t237) + t138 * t213 - t137 * t214) * m(4) + (t21 * t354 + t23 * t351 + t24 * t349 + t25 * t347) * t212 + ((-t17 * t355 - t18 * t352 - t19 * t350 - t20 * t348) * t212 + t414 * (t17 * t310 - t187 * (-t127 * t77 + t170 * t94 - t73 * t333)) + t415 * (t18 * t309 - t194 * (-t128 * t78 + t170 * t95 - t74 * t331)) + t416 * (t19 * t308 - t197 * (-t129 * t79 + t170 * t96 - t75 * t329)) + t417 * (t20 * t307 - t200 * (-t130 * t80 + t170 * t97 - t76 * t327))) * t253;];
tauX  = t22;
