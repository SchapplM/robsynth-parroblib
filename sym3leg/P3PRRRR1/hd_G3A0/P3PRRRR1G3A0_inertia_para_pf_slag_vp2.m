% Calculate inertia matrix for parallel robot
% P3PRRRR1G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:21
% EndTime: 2020-03-09 21:02:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (615->124), mult. (792->250), div. (441->10), fcn. (1089->18), ass. (0->121)
t336 = 2 * Ifges(3,4);
t277 = cos(qJ(3,3));
t260 = 0.1e1 / t277;
t335 = t260 ^ 2;
t279 = cos(qJ(3,2));
t262 = 0.1e1 / t279;
t334 = t262 ^ 2;
t281 = cos(qJ(3,1));
t264 = 0.1e1 / t281;
t333 = t264 ^ 2;
t332 = Ifges(3,1) + Ifges(2,3);
t268 = legFrame(3,2);
t250 = sin(t268);
t253 = cos(t268);
t272 = sin(qJ(2,3));
t278 = cos(qJ(2,3));
t238 = -t250 * t278 + t253 * t272;
t257 = 0.1e1 / t272;
t331 = t238 * t257;
t239 = t272 * t250 + t253 * t278;
t330 = t239 * t257;
t269 = legFrame(2,2);
t251 = sin(t269);
t254 = cos(t269);
t274 = sin(qJ(2,2));
t280 = cos(qJ(2,2));
t240 = -t251 * t280 + t254 * t274;
t258 = 0.1e1 / t274;
t329 = t240 * t258;
t241 = t274 * t251 + t254 * t280;
t328 = t241 * t258;
t270 = legFrame(1,2);
t252 = sin(t270);
t255 = cos(t270);
t276 = sin(qJ(2,1));
t282 = cos(qJ(2,1));
t242 = -t252 * t282 + t255 * t276;
t259 = 0.1e1 / t276;
t327 = t242 * t259;
t243 = t276 * t252 + t255 * t282;
t326 = t243 * t259;
t271 = sin(qJ(3,3));
t247 = mrSges(3,1) * t271 + mrSges(3,2) * t277;
t325 = t247 * t260;
t273 = sin(qJ(3,2));
t248 = mrSges(3,1) * t273 + mrSges(3,2) * t279;
t324 = t248 * t262;
t275 = sin(qJ(3,1));
t249 = mrSges(3,1) * t275 + mrSges(3,2) * t281;
t323 = t249 * t264;
t322 = t257 * t260;
t321 = t257 * t271;
t320 = t258 * t262;
t319 = t258 * t273;
t318 = t259 * t264;
t317 = t259 * t275;
t283 = 0.1e1 / pkin(2);
t316 = t260 * t283;
t284 = t277 ^ 2;
t315 = 0.1e1 / t284 * t278;
t314 = t262 * t283;
t285 = t279 ^ 2;
t313 = 0.1e1 / t285 * t280;
t312 = t264 * t283;
t286 = t281 ^ 2;
t311 = 0.1e1 / t286 * t282;
t244 = Ifges(3,5) * t271 + Ifges(3,6) * t277;
t310 = t244 * t257 * t335;
t245 = Ifges(3,5) * t273 + Ifges(3,6) * t279;
t309 = t245 * t258 * t334;
t246 = Ifges(3,5) * t275 + Ifges(3,6) * t281;
t308 = t246 * t259 * t333;
t307 = t250 * t322;
t306 = t250 * t316;
t305 = t251 * t320;
t304 = t251 * t314;
t303 = t252 * t318;
t302 = t252 * t312;
t301 = t253 * t322;
t300 = t253 * t316;
t299 = t254 * t320;
t298 = t254 * t314;
t297 = t255 * t318;
t296 = t255 * t312;
t295 = t260 * t321;
t294 = t262 * t319;
t293 = t264 * t317;
t292 = t283 * t315;
t291 = t283 * t313;
t290 = t283 * t311;
t289 = t315 * t321;
t288 = t313 * t319;
t287 = t311 * t317;
t267 = mrSges(2,2) - mrSges(3,3);
t266 = -Ifges(3,1) + Ifges(3,2);
t256 = m(1) + m(2) + m(3);
t237 = t275 * t281 * t336 + t266 * t286 + t332;
t236 = t273 * t279 * t336 + t266 * t285 + t332;
t235 = t271 * t277 * t336 + t266 * t284 + t332;
t234 = (t281 * mrSges(3,1) - t275 * mrSges(3,2) + mrSges(2,1)) * t282 - t276 * t267;
t233 = (t279 * mrSges(3,1) - t273 * mrSges(3,2) + mrSges(2,1)) * t280 - t274 * t267;
t232 = (t277 * mrSges(3,1) - t271 * mrSges(3,2) + mrSges(2,1)) * t278 - t272 * t267;
t231 = (-t234 * t296 + t243 * t256) * t259;
t230 = (t234 * t302 + t242 * t256) * t259;
t229 = (-t233 * t298 + t241 * t256) * t258;
t228 = (t233 * t304 + t240 * t256) * t258;
t227 = (-t232 * t300 + t239 * t256) * t257;
t226 = (t232 * t306 + t238 * t256) * t257;
t225 = -t276 * t249 * t312 + (-t234 * t290 + t256 * t264) * t317;
t224 = -t274 * t248 * t314 + (-t233 * t291 + t256 * t262) * t319;
t223 = -t272 * t247 * t316 + (-t232 * t292 + t256 * t260) * t321;
t222 = (t234 * t243 - t237 * t296) * t259;
t221 = (t234 * t242 + t237 * t302) * t259;
t220 = (t233 * t241 - t236 * t298) * t258;
t219 = (t233 * t240 + t236 * t304) * t258;
t218 = (t232 * t239 - t235 * t300) * t257;
t217 = (t232 * t238 + t235 * t306) * t257;
t216 = t246 * t312 + (t234 * t264 - t237 * t290) * t317;
t215 = t245 * t314 + (t233 * t262 - t236 * t291) * t319;
t214 = t244 * t316 + (t232 * t260 - t235 * t292) * t321;
t1 = [t227 * t330 + t229 * t328 + t231 * t326 + m(4) + (-t218 * t301 - t220 * t299 - t222 * t297) * t283, t227 * t331 + t229 * t329 + t231 * t327 + (t218 * t307 + t220 * t305 + t222 * t303) * t283, t227 * t295 + t229 * t294 + t231 * t293 + (-t218 * t289 - t220 * t288 - t222 * t287 - t239 * t325 - t241 * t324 - t243 * t323 + (-t253 * t310 - t254 * t309 - t255 * t308) * t283) * t283; t226 * t330 + t228 * t328 + t230 * t326 + (-t217 * t301 - t219 * t299 - t221 * t297) * t283, t226 * t331 + t228 * t329 + t230 * t327 + m(4) + (t217 * t307 + t219 * t305 + t221 * t303) * t283, t226 * t295 + t228 * t294 + t230 * t293 + (-t217 * t289 - t219 * t288 - t221 * t287 - t238 * t325 - t240 * t324 - t242 * t323 + (t250 * t310 + t251 * t309 + t252 * t308) * t283) * t283; t223 * t330 + t224 * t328 + t225 * t326 + (-t214 * t301 - t215 * t299 - t216 * t297) * t283, t223 * t331 + t224 * t329 + t225 * t327 + (t214 * t307 + t215 * t305 + t216 * t303) * t283, t223 * t295 + t224 * t294 + t225 * t293 + m(4) + (-t214 * t289 - t215 * t288 - t216 * t287 - t271 * t335 * t247 - t273 * t334 * t248 - t275 * t333 * t249 + ((Ifges(3,3) * t264 - t246 * t287) * t264 + (Ifges(3,3) * t262 - t245 * t288) * t262 + (Ifges(3,3) * t260 - t244 * t289) * t260) * t283) * t283;];
MX  = t1;
