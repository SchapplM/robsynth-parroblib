% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:05:53
% EndTime: 2022-11-07 13:06:00
% DurationCPUTime: 8.01s
% Computational Cost: add. (31869->598), mult. (55041->930), div. (4491->12), fcn. (38313->50), ass. (0->370)
t214 = (qJ(3,1) + pkin(5));
t140 = mrSges(3,2) * t214 - Ifges(3,6);
t209 = sin(pkin(7));
t210 = cos(pkin(7));
t398 = mrSges(3,1) * t214 - Ifges(3,5);
t415 = -t140 * t209 + t398 * t210;
t213 = (qJ(3,2) + pkin(5));
t139 = mrSges(3,2) * t213 - Ifges(3,6);
t399 = mrSges(3,1) * t213 - Ifges(3,5);
t414 = -t139 * t209 + t399 * t210;
t212 = (qJ(3,3) + pkin(5));
t138 = mrSges(3,2) * t212 - Ifges(3,6);
t400 = mrSges(3,1) * t212 - Ifges(3,5);
t413 = -t138 * t209 + t400 * t210;
t197 = t210 ^ 2;
t412 = 0.2e1 * t197;
t397 = 2 * pkin(1);
t317 = 2 * pkin(3);
t411 = pkin(2) / 0.2e1;
t410 = t397 / 0.2e1;
t203 = qJ(2,1) + pkin(7);
t168 = cos(t203);
t151 = pkin(3) * t168;
t234 = cos(qJ(2,1));
t189 = t234 * pkin(2);
t319 = t189 + t151;
t409 = pkin(1) + t319;
t202 = qJ(2,2) + pkin(7);
t167 = cos(t202);
t150 = pkin(3) * t167;
t232 = cos(qJ(2,2));
t188 = t232 * pkin(2);
t320 = t188 + t150;
t408 = pkin(1) + t320;
t201 = qJ(2,3) + pkin(7);
t166 = cos(t201);
t149 = pkin(3) * t166;
t230 = cos(qJ(2,3));
t187 = t230 * pkin(2);
t321 = t187 + t149;
t407 = pkin(1) + t321;
t169 = t209 * mrSges(3,1);
t215 = Ifges(3,2) - Ifges(3,1);
t101 = (-pkin(2) * mrSges(3,2) - t209 * t215) * t210 - pkin(2) * t169 + Ifges(3,4) * t412 + Ifges(2,4) - Ifges(3,4);
t171 = (t212 * m(3));
t217 = (mrSges(2,3) + mrSges(3,3));
t289 = m(2) * pkin(5) - mrSges(1,2) + t217;
t403 = -t289 - t171;
t172 = (t213 * m(3));
t402 = -t289 - t172;
t173 = (t214 * m(3));
t401 = -t289 - t173;
t122 = 0.1e1 / t321;
t218 = legFrame(3,2);
t177 = sin(t218);
t180 = cos(t218);
t239 = xDP(2);
t240 = xDP(1);
t85 = (t177 * t240 + t180 * t239) * t122;
t82 = t85 ^ 2;
t123 = 0.1e1 / t320;
t219 = legFrame(2,2);
t178 = sin(t219);
t181 = cos(t219);
t86 = (t178 * t240 + t181 * t239) * t123;
t83 = t86 ^ 2;
t124 = 0.1e1 / t319;
t220 = legFrame(1,2);
t179 = sin(t220);
t182 = cos(t220);
t87 = (t179 * t240 + t182 * t239) * t124;
t84 = t87 ^ 2;
t396 = 0.2e1 * mrSges(3,1);
t395 = 2 * mrSges(3,3);
t259 = pkin(3) ^ 2;
t367 = pkin(3) * t210;
t161 = pkin(2) * t367;
t260 = pkin(2) ^ 2;
t211 = t260 / 0.2e1;
t322 = t161 + t211;
t394 = -0.4e1 * pkin(1) * (t259 / 0.2e1 + t322);
t393 = 0.2e1 * t161;
t205 = t230 ^ 2;
t392 = 0.2e1 * t205;
t206 = t232 ^ 2;
t391 = 0.2e1 * t206;
t207 = t234 ^ 2;
t390 = 0.2e1 * t207;
t389 = 0.2e1 * t210;
t224 = sin(qJ(2,3));
t388 = -0.2e1 * t224;
t226 = sin(qJ(2,2));
t387 = -0.2e1 * t226;
t228 = sin(qJ(2,1));
t386 = -0.2e1 * t228;
t385 = -0.4e1 * t230;
t384 = -0.4e1 * t232;
t383 = -0.4e1 * t234;
t382 = -4 * pkin(5) - 4 * pkin(6);
t246 = m(3) * pkin(2);
t381 = mrSges(2,2) * pkin(1);
t154 = pkin(2) + t367;
t328 = t209 * t224;
t276 = pkin(3) * t328 - t154 * t230;
t108 = 0.1e1 / t276;
t194 = pkin(6) + t212;
t174 = 1 / t194;
t231 = cos(qJ(1,3));
t238 = xDP(3);
t225 = sin(qJ(1,3));
t368 = pkin(3) * t209;
t305 = t225 * t368;
t308 = t180 * t368;
t334 = t177 * t225;
t337 = t154 * t180;
t76 = (-t154 * t334 + t308) * t230 + (t177 * t305 + t337) * t224;
t311 = t177 * t368;
t331 = t180 * t225;
t340 = t154 * t177;
t79 = (t154 * t331 + t311) * t230 + t224 * (-t180 * t305 + t340);
t64 = (t231 * t238 - (t239 * t76 + t240 * t79) * t108) * t174;
t62 = pkin(1) * t64;
t327 = t209 * t226;
t275 = pkin(3) * t327 - t154 * t232;
t109 = 0.1e1 / t275;
t195 = pkin(6) + t213;
t175 = 1 / t195;
t233 = cos(qJ(1,2));
t227 = sin(qJ(1,2));
t304 = t227 * t368;
t307 = t181 * t368;
t333 = t178 * t227;
t336 = t154 * t181;
t77 = (-t154 * t333 + t307) * t232 + (t178 * t304 + t336) * t226;
t310 = t178 * t368;
t330 = t181 * t227;
t339 = t154 * t178;
t80 = (t154 * t330 + t310) * t232 + t226 * (-t181 * t304 + t339);
t65 = (t233 * t238 - (t239 * t77 + t240 * t80) * t109) * t175;
t63 = pkin(1) * t65;
t287 = pkin(1) * t225 - t194 * t231;
t302 = pkin(3) * (t197 - 0.1e1);
t380 = pkin(3) * (t225 * t302 + t287 * t328);
t286 = pkin(1) * t227 - t195 * t233;
t379 = pkin(3) * (t227 * t302 + t286 * t327);
t229 = sin(qJ(1,1));
t196 = pkin(6) + t214;
t235 = cos(qJ(1,1));
t285 = pkin(1) * t229 - t196 * t235;
t326 = t209 * t228;
t378 = pkin(3) * (t229 * t302 + t285 * t326);
t377 = pkin(5) * mrSges(2,2);
t274 = pkin(3) * t326 - t154 * t234;
t110 = 0.1e1 / t274;
t176 = 1 / t196;
t303 = t229 * t368;
t306 = t182 * t368;
t332 = t179 * t229;
t335 = t154 * t182;
t78 = (-t154 * t332 + t306) * t234 + (t179 * t303 + t335) * t228;
t309 = t179 * t368;
t329 = t182 * t229;
t338 = t154 * t179;
t81 = (t154 * t329 + t309) * t234 + t228 * (-t182 * t303 + t338);
t66 = (t235 * t238 - (t239 * t78 + t240 * t81) * t110) * t176;
t61 = t66 * pkin(1);
t376 = t122 / 0.2e1;
t375 = t123 / 0.2e1;
t374 = t124 / 0.2e1;
t373 = m(3) * t260;
t183 = mrSges(2,1) + t246;
t372 = pkin(1) * t183;
t371 = pkin(1) * t224;
t370 = pkin(1) * t226;
t369 = pkin(1) * t228;
t366 = pkin(3) * t260;
t365 = t259 * pkin(2);
t364 = t64 * t85;
t363 = t65 * t86;
t362 = t66 * t87;
t170 = mrSges(3,1) * t210;
t361 = mrSges(3,2) * t209;
t360 = mrSges(3,2) * t210;
t358 = Ifges(3,4) * t209;
t325 = t215 * t197;
t91 = 0.2e1 * t325 + (pkin(2) * t396 + 0.4e1 * t358) * t210 + t373 - 0.2e1 * pkin(2) * t361 + Ifges(2,2) - Ifges(2,1) - t215;
t357 = t224 * t91;
t356 = t226 * t91;
t355 = t228 * t91;
t354 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t353 = -t183 * pkin(5) + Ifges(2,5);
t148 = t169 + mrSges(2,2);
t352 = t101 * t205;
t351 = t101 * t206;
t350 = t101 * t207;
t349 = t108 * t174;
t348 = t109 * t175;
t347 = t110 * t176;
t346 = t122 * t177;
t345 = t122 * t180;
t344 = t123 * t178;
t343 = t123 * t181;
t342 = t124 * t179;
t341 = t124 * t182;
t324 = -0.2e1 * t246;
t323 = pkin(2) * t317;
t318 = t259 + t260;
t316 = -0.2e1 * t366;
t315 = -0.2e1 * t148 * pkin(1);
t314 = -0.2e1 * t365;
t105 = t305 * t388 + t287;
t112 = (t197 - 0.1e1 / 0.2e1) * t259 + t322;
t132 = -t368 + t371;
t162 = pkin(1) * t368;
t268 = t259 * t412 - t259 + t260 + t393;
t98 = t224 * t268 + t162;
t55 = (-t112 * t334 + t154 * t308) * t392 + (-t105 * t340 + t180 * t98) * t230 + t177 * t380 + t132 * t337;
t40 = t55 * t239 * t349;
t56 = (t112 * t331 + t154 * t311) * t392 + (t105 * t337 + t177 * t98) * t230 - t180 * t380 + t132 * t340;
t41 = t56 * t240 * t349;
t102 = t225 * t194 + t407 * t231;
t92 = t102 * t174 * t238;
t28 = -t41 - t40 + t92;
t106 = t304 * t387 + t286;
t133 = -t368 + t370;
t99 = t226 * t268 + t162;
t57 = (-t112 * t333 + t154 * t307) * t391 + (-t106 * t339 + t181 * t99) * t232 + t178 * t379 + t133 * t336;
t42 = t57 * t239 * t348;
t58 = (t112 * t330 + t154 * t310) * t391 + (t106 * t336 + t178 * t99) * t232 - t181 * t379 + t133 * t339;
t43 = t58 * t240 * t348;
t103 = t227 * t195 + t408 * t233;
t93 = t103 * t175 * t238;
t29 = -t43 - t42 + t93;
t100 = t228 * t268 + t162;
t107 = t303 * t386 + t285;
t134 = -t368 + t369;
t59 = (-t112 * t332 + t154 * t306) * t390 + (t100 * t182 - t107 * t338) * t234 + t179 * t378 + t134 * t335;
t44 = t59 * t239 * t347;
t60 = (t112 * t329 + t154 * t309) * t390 + (t100 * t179 + t107 * t335) * t234 - t182 * t378 + t134 * t338;
t45 = t60 * t240 * t347;
t104 = t229 * t196 + t409 * t235;
t94 = t104 * t176 * t238;
t30 = -t45 - t44 + t94;
t292 = -m(3) * qJ(3,3) - mrSges(3,3);
t295 = Ifges(2,6) - t377;
t73 = (pkin(2) * t292 + t353 - t413) * t224 + (-t138 * t210 - t209 * t400 + t295) * t230;
t301 = t73 * t349;
t293 = -m(3) * qJ(3,2) - mrSges(3,3);
t74 = (pkin(2) * t293 + t353 - t414) * t226 + (-t139 * t210 - t209 * t399 + t295) * t232;
t300 = t74 * t348;
t294 = -m(3) * qJ(3,1) - mrSges(3,3);
t75 = (pkin(2) * t294 + t353 - t415) * t228 + (-t140 * t210 - t209 * t398 + t295) * t234;
t299 = t75 * t347;
t126 = pkin(2) * t224 + pkin(3) * sin(t201);
t298 = t122 * t126 * t82;
t128 = pkin(2) * t226 + pkin(3) * sin(t202);
t297 = t123 * t128 * t83;
t127 = pkin(2) * t228 + pkin(3) * sin(t203);
t296 = t124 * t127 * t84;
t291 = -(2 * pkin(3) * t259) - 0.4e1 * t366;
t290 = t246 - t361;
t288 = pkin(1) * t396 * t210 + (mrSges(2,1) + t290) * t397;
t114 = -t170 - t290;
t284 = -t361 + t170;
t283 = t169 + t360;
t118 = g(1) * t180 - g(2) * t177;
t282 = -g(3) * t231 - t118 * t225;
t119 = g(1) * t181 - g(2) * t178;
t281 = -g(3) * t233 - t119 * t227;
t120 = g(1) * t182 - g(2) * t179;
t280 = -g(3) * t235 - t120 * t229;
t113 = -mrSges(2,1) + t114;
t129 = t148 + t360;
t279 = -t113 * t230 - t129 * t224;
t278 = -t113 * t232 - t129 * t226;
t277 = -t113 * t234 - t129 * t228;
t257 = pkin(5) ^ 2;
t261 = pkin(1) ^ 2;
t273 = -(2 * t257) - (2 * t261) - t318 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t272 = 0.2e1 * t101;
t271 = t283 * t230;
t270 = t283 * t232;
t269 = t283 * t234;
t266 = (mrSges(2,2) + t283) * t397;
t264 = (m(2) * t257) + (2 * t217 * pkin(5)) + (m(2) + m(3)) * t261 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t325;
t256 = 0.2e1 * qJ(2,1);
t255 = 0.2e1 * qJ(2,2);
t254 = 0.2e1 * qJ(2,3);
t253 = -pkin(5) / 0.2e1;
t252 = 0.2e1 * pkin(7);
t247 = m(3) * pkin(1);
t245 = m(3) * pkin(5);
t244 = Ifges(3,5) / 0.2e1;
t243 = Ifges(3,6) / 0.2e1;
t223 = xDDP(1);
t222 = xDDP(2);
t221 = xDDP(3);
t200 = t254 + pkin(7);
t199 = pkin(7) + t256;
t198 = pkin(7) + t255;
t186 = -qJ(3,1) / 0.2e1 + t253;
t185 = -qJ(3,2) / 0.2e1 + t253;
t184 = -qJ(3,3) / 0.2e1 + t253;
t160 = pkin(2) * t211 + t365;
t158 = -t377 / 0.2e1 + Ifges(2,6) / 0.2e1;
t157 = 0.2e1 * t203;
t156 = 0.2e1 * t202;
t155 = 0.2e1 * t201;
t147 = (m(2) * pkin(1)) + mrSges(1,1) + t247;
t146 = t245 - t294;
t145 = t245 - t293;
t144 = t245 - t292;
t137 = t147 * g(3);
t135 = t393 + t318;
t117 = g(1) * t179 + g(2) * t182;
t116 = g(1) * t178 + g(2) * t181;
t115 = g(1) * t177 + g(2) * t180;
t111 = 0.2e1 * pkin(2) * t284 + Ifges(2,3) + Ifges(3,3) + t373;
t90 = t114 * t234 + t228 * t283 - t247;
t89 = t114 * t232 + t226 * t283 - t247;
t88 = t114 * t230 + t224 * t283 - t247;
t72 = (m(3) * t104 + t235 * t90) * t176;
t71 = (m(3) * t103 + t233 * t89) * t175;
t70 = (m(3) * t102 + t231 * t88) * t174;
t69 = t91 * t207 + (t228 * t272 + t288) * t234 + (-mrSges(3,2) * t369 - t358) * t389 + t228 * t315 + (t214 ^ 2) * m(3) + qJ(3,1) * t395 + t264;
t68 = t91 * t206 + (t226 * t272 + t288) * t232 + (-mrSges(3,2) * t370 - t358) * t389 + t226 * t315 + (t213 ^ 2) * m(3) + qJ(3,2) * t395 + t264;
t67 = t91 * t205 + (t224 * t272 + t288) * t230 + (-mrSges(3,2) * t371 - t358) * t389 + t224 * t315 + (t212 ^ 2) * m(3) + qJ(3,3) * t395 + t264;
t54 = t66 * t372;
t53 = t65 * t372;
t52 = t64 * t372;
t51 = t111 * t342 - t299 * t81;
t50 = t111 * t344 - t300 * t80;
t49 = t111 * t346 - t301 * t79;
t48 = t111 * t341 - t299 * t78;
t47 = t111 * t343 - t300 * t77;
t46 = t111 * t345 - t301 * t76;
t39 = (t104 * t90 + t235 * t69) * t176;
t38 = (t103 * t89 + t233 * t68) * t175;
t37 = (t102 * t88 + t231 * t67) * t174;
t36 = (m(3) * t60 + t81 * t90) * t347;
t35 = (m(3) * t58 + t80 * t89) * t348;
t34 = (m(3) * t56 + t79 * t88) * t349;
t33 = (m(3) * t59 + t78 * t90) * t347;
t32 = (m(3) * t57 + t77 * t89) * t348;
t31 = (m(3) * t55 + t76 * t88) * t349;
t27 = t75 * t342 - (t60 * t90 + t69 * t81) * t347;
t26 = t74 * t344 - (t58 * t89 + t68 * t80) * t348;
t25 = t73 * t346 - (t56 * t88 + t67 * t79) * t349;
t24 = t75 * t341 - (t59 * t90 + t69 * t78) * t347;
t23 = t74 * t343 - (t57 * t89 + t68 * t77) * t348;
t22 = t73 * t345 - (t55 * t88 + t67 * t76) * t349;
t18 = t63 + t43 / 0.2e1 + t42 / 0.2e1 - t93 / 0.2e1;
t17 = t62 + t41 / 0.2e1 + t40 / 0.2e1 - t92 / 0.2e1;
t16 = t61 + t45 / 0.2e1 + t44 / 0.2e1 - t94 / 0.2e1;
t15 = (-t84 * t135 / (t189 + (t210 * t234 - t326) * pkin(3)) + (t274 * t66 + 0.2e1 * t30 - t61) * t66) * t176;
t14 = (-t83 * t135 / (t188 + (t210 * t232 - t327) * pkin(3)) + (t275 * t65 + 0.2e1 * t29 - t63) * t65) * t175;
t13 = (-t82 * t135 / (t187 + (t210 * t230 - t328) * pkin(3)) + (t276 * t64 + 0.2e1 * t28 - t62) * t64) * t174;
t12 = ((cos(qJ(2,1) - pkin(7)) * t316 + cos(qJ(2,1) + t252) * t314 + t291 * t168 + t160 * t383 + t394) * t84 * t374 + (-0.2e1 * t16 * t151 + (-t260 * cos(t256) - t259 * cos(t157) + t273 + (t382 - 0.2e1 * qJ(3,1)) * qJ(3,1)) * t66 / 0.2e1 + (t16 * t383 + (-cos(t199) - t210) * t66 * t317) * t411 + (t127 + (sin(t199) * t323 + sin(t157) * t259 + sin(t256) * t260) * t374) * t87 * t196 + (t409 + t410) * t30) * t66) * t176;
t11 = ((cos(qJ(2,2) - pkin(7)) * t316 + cos(t252 + qJ(2,2)) * t314 + t291 * t167 + t160 * t384 + t394) * t83 * t375 + (-0.2e1 * t18 * t150 + (-t259 * cos(t156) - t260 * cos(t255) + t273 + (t382 - 0.2e1 * qJ(3,2)) * qJ(3,2)) * t65 / 0.2e1 + (t18 * t384 + (-cos(t198) - t210) * t65 * t317) * t411 + (t128 + (sin(t198) * t323 + sin(t156) * t259 + sin(t255) * t260) * t375) * t86 * t195 + (t408 + t410) * t29) * t65) * t175;
t10 = ((cos(qJ(2,3) - pkin(7)) * t316 + cos(t252 + qJ(2,3)) * t314 + t291 * t166 + t160 * t385 + t394) * t82 * t376 + (-0.2e1 * t17 * t149 + (-t259 * cos(t155) - t260 * cos(t254) + t273 + (t382 - 0.2e1 * qJ(3,3)) * qJ(3,3)) * t64 / 0.2e1 + (t17 * t385 + (-cos(t200) - t210) * t64 * t317) * t411 + (t126 + (sin(t200) * t323 + sin(t155) * t259 + sin(t254) * t260) * t376) * t85 * t194 + (t407 + t410) * t28) * t64) * t174;
t9 = -t75 * t15 + t111 * t296 + t66 * ((t30 * t324 + t54) * t228 + (t228 * t284 + t269) * (t61 + 0.2e1 * t45 + 0.2e1 * t44 - 0.2e1 * t94) + (-0.2e1 * t350 + (t355 + t381) * t234 + t101) * t66) + (t113 * t280 + t117 * t129) * t228 - (-t117 * t113 + t129 * t280) * t234;
t8 = -t74 * t14 + t111 * t297 + t65 * ((t29 * t324 + t53) * t226 + (t226 * t284 + t270) * (t63 + 0.2e1 * t43 + 0.2e1 * t42 - 0.2e1 * t93) + (-0.2e1 * t351 + (t356 + t381) * t232 + t101) * t65) + (t113 * t281 + t116 * t129) * t226 - (-t116 * t113 + t129 * t281) * t232;
t7 = -t73 * t13 + t111 * t298 + t64 * ((t28 * t324 + t52) * t224 + (t224 * t284 + t271) * (t62 + 0.2e1 * t41 + 0.2e1 * t40 - 0.2e1 * t92) + (-0.2e1 * t352 + (t357 + t381) * t230 + t101) * t64) + (t113 * t282 + t115 * t129) * t224 - (-t115 * t113 + t129 * t282) * t230;
t6 = -t90 * t15 - (t173 + mrSges(3,3)) * t66 ^ 2 + (-g(3) * t229 + t120 * t235 - t12) * m(3) + 0.2e1 * (-t114 * t228 + t269) * t362;
t5 = -t89 * t14 - (t172 + mrSges(3,3)) * t65 ^ 2 + (-g(3) * t227 + t119 * t233 - t11) * m(3) + 0.2e1 * (-t114 * t226 + t270) * t363;
t4 = -t88 * t13 - (t171 + mrSges(3,3)) * t64 ^ 2 + (-g(3) * t225 + t118 * t231 - t10) * m(3) + 0.2e1 * (-t114 * t224 + t271) * t364;
t3 = -t69 * t15 + t75 * t296 - t90 * t12 + 0.4e1 * t350 * t362 - ((t146 * pkin(2) + t354 + t415) * t87 + (t266 + 0.2e1 * t355) * t66) * t87 * t234 + (((mrSges(3,2) * t186 + t243) * t87 + mrSges(3,1) * t61) * t210 + ((mrSges(3,1) * t186 + t244) * t87 - mrSges(3,2) * t61) * t209 + t158 * t87 + t54) * t87 * t386 + 0.2e1 * (-t101 * t87 + t30 * t146) * t66 + (t401 * g(3) + (-t147 - t277) * t120) * t235 + (t277 * g(3) + t401 * t120 + t137) * t229;
t2 = -t68 * t14 + t74 * t297 - t89 * t11 + 0.4e1 * t351 * t363 - ((t145 * pkin(2) + t354 + t414) * t86 + (t266 + 0.2e1 * t356) * t65) * t86 * t232 + (((mrSges(3,2) * t185 + t243) * t86 + mrSges(3,1) * t63) * t210 + ((mrSges(3,1) * t185 + t244) * t86 - mrSges(3,2) * t63) * t209 + t158 * t86 + t53) * t86 * t387 + 0.2e1 * (-t101 * t86 + t29 * t145) * t65 + (t402 * g(3) + (-t147 - t278) * t119) * t233 + (t278 * g(3) + t402 * t119 + t137) * t227;
t1 = -t67 * t13 + t73 * t298 - t88 * t10 + 0.4e1 * t352 * t364 - ((t144 * pkin(2) + t354 + t413) * t85 + (t266 + 0.2e1 * t357) * t64) * t85 * t230 + (((mrSges(3,2) * t184 + t243) * t85 + mrSges(3,1) * t62) * t210 + ((mrSges(3,1) * t184 + t244) * t85 - mrSges(3,2) * t62) * t209 + t158 * t85 + t52) * t85 * t388 + 0.2e1 * (-t101 * t85 + t28 * t144) * t64 + (t403 * g(3) + (-t147 - t279) * t118) * t231 + (t279 * g(3) + t403 * t118 + t137) * t225;
t19 = [t7 * t346 + t8 * t344 + t9 * t342 - g(1) * m(4) + (t341 * t51 + t343 * t50 + t345 * t49) * t222 + (t342 * t51 + t344 * t50 + t346 * t49 + m(4)) * t223 + ((-t104 * t36 + t235 * t27) * t221 - ((t27 * t81 - t36 * t60) * t223 + (t27 * t78 - t36 * t59) * t222 + t81 * t3 + t60 * t6) * t110) * t176 + ((-t103 * t35 + t233 * t26) * t221 - ((t26 * t80 - t35 * t58) * t223 + (t26 * t77 - t35 * t57) * t222 + t80 * t2 + t58 * t5) * t109) * t175 + ((-t102 * t34 + t231 * t25) * t221 - ((t25 * t79 - t34 * t56) * t223 + (t25 * t76 - t34 * t55) * t222 + t79 * t1 + t56 * t4) * t108) * t174; t7 * t345 + t8 * t343 + t9 * t341 - g(2) * m(4) + (t342 * t48 + t344 * t47 + t346 * t46) * t223 + (t341 * t48 + t343 * t47 + t345 * t46 + m(4)) * t222 + ((-t104 * t33 + t235 * t24) * t221 - ((t24 * t81 - t33 * t60) * t223 + (t24 * t78 - t33 * t59) * t222 + t78 * t3 + t59 * t6) * t110) * t176 + ((-t103 * t32 + t23 * t233) * t221 - ((t23 * t80 - t32 * t58) * t223 + (t23 * t77 - t32 * t57) * t222 + t77 * t2 + t57 * t5) * t109) * t175 + ((-t102 * t31 + t22 * t231) * t221 - ((t22 * t79 - t31 * t56) * t223 + (t22 * t76 - t31 * t55) * t222 + t76 * t1 + t55 * t4) * t108) * t174; (-g(3) + t221) * m(4) + ((t72 * t221 + t6) * t104 + (t39 * t221 + t3 + (t179 * t223 + t182 * t222) * t75 * t124) * t235 - ((t39 * t81 + t60 * t72) * t223 + (t39 * t78 + t59 * t72) * t222) * t110) * t176 + ((t71 * t221 + t5) * t103 + (t38 * t221 + t2 + (t178 * t223 + t181 * t222) * t74 * t123) * t233 - ((t38 * t80 + t58 * t71) * t223 + (t38 * t77 + t57 * t71) * t222) * t109) * t175 + ((t70 * t221 + t4) * t102 + (t37 * t221 + t1 + (t177 * t223 + t180 * t222) * t73 * t122) * t231 - ((t37 * t79 + t56 * t70) * t223 + (t37 * t76 + t55 * t70) * t222) * t108) * t174;];
tauX  = t19;
