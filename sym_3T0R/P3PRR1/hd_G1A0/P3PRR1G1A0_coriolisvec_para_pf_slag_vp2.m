% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRR1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:37
% EndTime: 2019-05-03 14:47:37
% DurationCPUTime: 0.60s
% Computational Cost: add. (640->122), mult. (1214->257), div. (462->4), fcn. (1318->14), ass. (0->115)
t286 = xDP(3);
t276 = t286 ^ 2;
t289 = m(1) + m(2);
t324 = pkin(2) * t289;
t290 = xP(3);
t271 = sin(t290);
t272 = cos(t290);
t293 = koppelP(3,2);
t296 = koppelP(3,1);
t250 = t271 * t296 + t272 * t293;
t253 = -t271 * t293 + t272 * t296;
t277 = legFrame(3,3);
t262 = sin(t277);
t265 = cos(t277);
t287 = xDP(2);
t288 = xDP(1);
t280 = sin(qJ(2,3));
t273 = 0.1e1 / t280;
t299 = 0.1e1 / pkin(2);
t312 = t273 * t299;
t226 = (-t265 * (-t250 * t286 + t288) - t262 * (t253 * t286 + t287)) * t312;
t223 = t226 ^ 2;
t323 = t223 * t273;
t294 = koppelP(2,2);
t297 = koppelP(2,1);
t251 = t271 * t297 + t272 * t294;
t254 = -t271 * t294 + t272 * t297;
t278 = legFrame(2,3);
t263 = sin(t278);
t266 = cos(t278);
t281 = sin(qJ(2,2));
t274 = 0.1e1 / t281;
t309 = t274 * t299;
t227 = (-t266 * (-t251 * t286 + t288) - t263 * (t254 * t286 + t287)) * t309;
t224 = t227 ^ 2;
t322 = t224 * t274;
t295 = koppelP(1,2);
t298 = koppelP(1,1);
t252 = t271 * t298 + t272 * t295;
t255 = -t271 * t295 + t272 * t298;
t279 = legFrame(1,3);
t264 = sin(t279);
t267 = cos(t279);
t282 = sin(qJ(2,1));
t275 = 0.1e1 / t282;
t306 = t275 * t299;
t228 = (-t267 * (-t252 * t286 + t288) - t264 * (t255 * t286 + t287)) * t306;
t225 = t228 ^ 2;
t321 = t225 * t275;
t320 = t262 * t299;
t319 = t263 * t299;
t318 = t264 * t299;
t317 = t265 * t299;
t316 = t266 * t299;
t315 = t267 * t299;
t283 = cos(qJ(2,3));
t259 = mrSges(2,1) * t283 - t280 * mrSges(2,2);
t209 = (-t259 * t283 + t324) * t323;
t314 = t273 * t209;
t313 = t273 * t276;
t284 = cos(qJ(2,2));
t260 = mrSges(2,1) * t284 - t281 * mrSges(2,2);
t210 = (-t260 * t284 + t324) * t322;
t311 = t274 * t210;
t310 = t274 * t276;
t285 = cos(qJ(2,1));
t261 = mrSges(2,1) * t285 - t282 * mrSges(2,2);
t208 = (-t261 * t285 + t324) * t321;
t308 = t275 * t208;
t307 = t275 * t276;
t256 = mrSges(2,1) * t280 + mrSges(2,2) * t283;
t305 = t256 * t323;
t257 = mrSges(2,1) * t281 + mrSges(2,2) * t284;
t304 = t257 * t322;
t258 = mrSges(2,1) * t282 + mrSges(2,2) * t285;
t303 = t258 * t321;
t211 = (-Ifges(2,3) * t283 + pkin(2) * t259) * t323;
t302 = t211 * t312;
t212 = (-Ifges(2,3) * t284 + pkin(2) * t260) * t322;
t301 = t212 * t309;
t213 = (-Ifges(2,3) * t285 + pkin(2) * t261) * t321;
t300 = t213 * t306;
t292 = mrSges(3,1);
t291 = mrSges(3,2);
t249 = t264 * t285 + t267 * t282;
t248 = -t264 * t282 + t267 * t285;
t247 = t263 * t284 + t266 * t281;
t246 = -t263 * t281 + t266 * t284;
t245 = t262 * t283 + t265 * t280;
t244 = -t262 * t280 + t265 * t283;
t243 = (-Ifges(2,3) * t318 + t249 * t261) * t275;
t242 = (-Ifges(2,3) * t315 + t248 * t261) * t275;
t241 = (-Ifges(2,3) * t319 + t247 * t260) * t274;
t240 = (-Ifges(2,3) * t316 + t246 * t260) * t274;
t239 = (-Ifges(2,3) * t320 + t245 * t259) * t273;
t238 = (-Ifges(2,3) * t317 + t244 * t259) * t273;
t237 = (t249 * t289 - t261 * t318) * t275;
t236 = (t248 * t289 - t261 * t315) * t275;
t235 = (t247 * t289 - t260 * t319) * t274;
t234 = (t246 * t289 - t260 * t316) * t274;
t233 = (t245 * t289 - t259 * t320) * t273;
t232 = (t244 * t289 - t259 * t317) * t273;
t231 = (t252 * t267 - t255 * t264) * t306;
t230 = (t251 * t266 - t254 * t263) * t309;
t229 = (t250 * t265 - t253 * t262) * t312;
t222 = (-t248 * t252 + t249 * t255) * t275;
t221 = (-t246 * t251 + t247 * t254) * t274;
t220 = (-t244 * t250 + t245 * t253) * t273;
t219 = Ifges(2,3) * t231 + t222 * t261;
t218 = Ifges(2,3) * t230 + t221 * t260;
t217 = Ifges(2,3) * t229 + t220 * t259;
t216 = t222 * t289 + t231 * t261;
t215 = t221 * t289 + t230 * t260;
t214 = t220 * t289 + t229 * t259;
t1 = [(-(t236 * t248 - t242 * t315) * t255 - (t236 * t249 - t242 * t318) * t252) * t307 + t248 * t308 - t267 * t300 - t248 * t303 + (-(t234 * t246 - t240 * t316) * t254 - (t234 * t247 - t240 * t319) * t251) * t310 + t246 * t311 - t266 * t301 - t246 * t304 + (-(t232 * t244 - t238 * t317) * t253 - (t232 * t245 - t238 * t320) * t250) * t313 + t244 * t314 - t265 * t302 - t244 * t305 - t276 * (-t271 * t291 + t272 * t292); (-(t237 * t248 - t243 * t315) * t255 - (t237 * t249 - t243 * t318) * t252) * t307 + t249 * t308 - t264 * t300 - t249 * t303 + (-(t235 * t246 - t241 * t316) * t254 - (t235 * t247 - t241 * t319) * t251) * t310 + t247 * t311 - t263 * t301 - t247 * t304 + (-(t233 * t244 - t239 * t317) * t253 - (t233 * t245 - t239 * t320) * t250) * t313 + t245 * t314 - t262 * t302 - t245 * t305 - t276 * (t271 * t292 + t272 * t291); (-(t216 * t248 - t219 * t315) * t255 - (t216 * t249 - t219 * t318) * t252) * t307 + t231 * t213 + (-(t215 * t246 - t218 * t316) * t254 - (t215 * t247 - t218 * t319) * t251) * t310 + t230 * t212 + (-(t214 * t244 - t217 * t317) * t253 - (t214 * t245 - t217 * t320) * t250) * t313 + t229 * t211 + (-t225 * t258 + t208) * t222 + (-t224 * t257 + t210) * t221 + (-t223 * t256 + t209) * t220;];
taucX  = t1;
