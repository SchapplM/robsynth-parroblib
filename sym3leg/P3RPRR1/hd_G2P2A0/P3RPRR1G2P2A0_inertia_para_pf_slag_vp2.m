% Calculate inertia matrix for parallel robot
% P3RPRR1G2P2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G2P2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:55
% EndTime: 2020-03-09 21:24:55
% DurationCPUTime: 0.68s
% Computational Cost: add. (2637->178), mult. (2346->234), div. (243->4), fcn. (1422->62), ass. (0->130)
t331 = sin(qJ(3,3));
t348 = pkin(7) + qJ(3,3);
t280 = 0.1e1 / (pkin(1) * sin(t348) + t331 * pkin(2));
t360 = t280 / 0.2e1;
t332 = sin(qJ(3,2));
t349 = pkin(7) + qJ(3,2);
t281 = 0.1e1 / (pkin(1) * sin(t349) + t332 * pkin(2));
t358 = t281 / 0.2e1;
t333 = sin(qJ(3,1));
t350 = pkin(7) + qJ(3,1);
t282 = 0.1e1 / (pkin(1) * sin(t350) + t333 * pkin(2));
t356 = t282 / 0.2e1;
t341 = 0.1e1 / pkin(3);
t392 = t341 / 0.2e1;
t330 = legFrame(1,2);
t353 = qJ(1,1) + pkin(7);
t308 = t330 + t353;
t302 = qJ(3,1) + t308;
t309 = -t330 + t353;
t303 = qJ(3,1) + t309;
t277 = -cos(t303) + cos(t302);
t391 = t277 * t356;
t329 = legFrame(2,2);
t352 = qJ(1,2) + pkin(7);
t306 = t329 + t352;
t300 = qJ(3,2) + t306;
t307 = -t329 + t352;
t301 = qJ(3,2) + t307;
t276 = -cos(t301) + cos(t300);
t390 = t276 * t358;
t328 = legFrame(3,2);
t351 = qJ(1,3) + pkin(7);
t304 = t328 + t351;
t298 = qJ(3,3) + t304;
t305 = -t328 + t351;
t299 = qJ(3,3) + t305;
t275 = -cos(t299) + cos(t298);
t389 = t275 * t360;
t274 = sin(t302) + sin(t303);
t388 = t274 * t356;
t273 = sin(t300) + sin(t301);
t387 = t273 * t358;
t272 = sin(t298) + sin(t299);
t386 = t272 * t360;
t326 = sin(pkin(7));
t385 = pkin(1) * t326;
t384 = Ifges(3,3) * t341;
t383 = m(3) * pkin(2) + mrSges(2,1);
t313 = qJ(1,3) + t328;
t314 = qJ(1,3) - t328;
t257 = t272 * pkin(3) + (sin(t304) + sin(t305)) * pkin(2) + (sin(t313) + sin(t314)) * pkin(1);
t382 = t257 * t280;
t315 = qJ(1,2) + t329;
t316 = qJ(1,2) - t329;
t258 = t273 * pkin(3) + (sin(t306) + sin(t307)) * pkin(2) + (sin(t315) + sin(t316)) * pkin(1);
t381 = t258 * t281;
t317 = qJ(1,1) + t330;
t318 = qJ(1,1) - t330;
t259 = t274 * pkin(3) + (sin(t308) + sin(t309)) * pkin(2) + (sin(t317) + sin(t318)) * pkin(1);
t380 = t259 * t282;
t260 = -t275 * pkin(3) + (cos(t305) - cos(t304)) * pkin(2) + (cos(t314) - cos(t313)) * pkin(1);
t379 = t260 * t280;
t261 = -t276 * pkin(3) + (cos(t307) - cos(t306)) * pkin(2) + (cos(t316) - cos(t315)) * pkin(1);
t378 = t261 * t281;
t262 = -t277 * pkin(3) + (cos(t309) - cos(t308)) * pkin(2) + (cos(t318) - cos(t317)) * pkin(1);
t377 = t262 * t282;
t327 = cos(pkin(7));
t346 = pkin(1) * t327 + pkin(2);
t278 = -mrSges(3,1) * t385 - t346 * mrSges(3,2);
t279 = t346 * mrSges(3,1) - mrSges(3,2) * t385;
t334 = cos(qJ(3,3));
t266 = t278 * t331 + t279 * t334 + Ifges(3,3);
t376 = t266 * t341;
t335 = cos(qJ(3,2));
t267 = t278 * t332 + t279 * t335 + Ifges(3,3);
t375 = t267 * t341;
t336 = cos(qJ(3,1));
t268 = t278 * t333 + t279 * t336 + Ifges(3,3);
t374 = t268 * t341;
t310 = cos(qJ(1,3) + t348);
t269 = -pkin(3) * t310 - pkin(2) * cos(t351) - cos(qJ(1,3)) * pkin(1);
t373 = t269 * t280;
t372 = t269 * t341;
t311 = cos(qJ(1,2) + t349);
t270 = -pkin(3) * t311 - pkin(2) * cos(t352) - cos(qJ(1,2)) * pkin(1);
t371 = t270 * t281;
t370 = t270 * t341;
t312 = cos(qJ(1,1) + t350);
t271 = -pkin(3) * t312 - pkin(2) * cos(t353) - cos(qJ(1,1)) * pkin(1);
t369 = t271 * t282;
t368 = t271 * t341;
t361 = t280 * t310;
t359 = t281 * t311;
t357 = t282 * t312;
t355 = 0.2e1 * pkin(1);
t354 = 0.2e1 * pkin(2);
t319 = sin(t328);
t320 = sin(t329);
t321 = sin(t330);
t322 = cos(t328);
t323 = cos(t329);
t324 = cos(t330);
t338 = m(2) + m(3);
t347 = (t319 * t322 + t320 * t323 + t321 * t324) * t338;
t345 = m(3) * pkin(2) ^ 2 + t338 * pkin(1) ^ 2 + Ifges(1,3) + Ifges(2,3) + Ifges(3,3);
t344 = t334 * mrSges(3,1) - mrSges(3,2) * t331;
t343 = t335 * mrSges(3,1) - mrSges(3,2) * t332;
t342 = t336 * mrSges(3,1) - mrSges(3,2) * t333;
t265 = t342 * t354 + ((t342 + t383) * t327 - (mrSges(3,1) * t333 + t336 * mrSges(3,2) + mrSges(2,2)) * t326) * t355 + t345;
t264 = t343 * t354 + ((t343 + t383) * t327 - (mrSges(3,1) * t332 + t335 * mrSges(3,2) + mrSges(2,2)) * t326) * t355 + t345;
t263 = t344 * t354 + ((t344 + t383) * t327 - (mrSges(3,1) * t331 + t334 * mrSges(3,2) + mrSges(2,2)) * t326) * t355 + t345;
t256 = (Ifges(3,3) * t368 + t268 * t312) * t282;
t255 = (Ifges(3,3) * t370 + t267 * t311) * t281;
t254 = (Ifges(3,3) * t372 + t266 * t310) * t280;
t253 = (t265 * t312 + t268 * t368) * t282;
t252 = (t264 * t311 + t267 * t370) * t281;
t251 = (t263 * t310 + t266 * t372) * t280;
t250 = (t262 * t384 + t268 * t277) * t356;
t249 = (t261 * t384 + t267 * t276) * t358;
t248 = (t260 * t384 + t266 * t275) * t360;
t247 = (-t259 * t384 + t268 * t274) * t356;
t246 = (-t258 * t384 + t267 * t273) * t358;
t245 = (-t257 * t384 + t266 * t272) * t360;
t244 = (t262 * t374 + t265 * t277) * t356;
t243 = (t261 * t375 + t264 * t276) * t358;
t242 = (t260 * t376 + t263 * t275) * t360;
t241 = (-t259 * t374 + t265 * t274) * t356;
t240 = (-t258 * t375 + t264 * t273) * t358;
t239 = (-t257 * t376 + t263 * t272) * t360;
t1 = [m(4) + (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) * t338 + t239 * t386 + t240 * t387 + t241 * t388 + (-t245 * t382 - t246 * t381 - t247 * t380) * t392, t239 * t389 + t240 * t390 + t241 * t391 + (t245 * t379 + t246 * t378 + t247 * t377) * t392 + t347, t239 * t361 + t240 * t359 + t241 * t357 + (t245 * t373 + t246 * t371 + t247 * t369) * t341; t242 * t386 + t243 * t387 + t244 * t388 + (-t248 * t382 - t249 * t381 - t250 * t380) * t392 + t347, m(4) + (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) * t338 + t242 * t389 + t243 * t390 + t244 * t391 + (t248 * t379 + t249 * t378 + t250 * t377) * t392, t242 * t361 + t243 * t359 + t244 * t357 + (t248 * t373 + t249 * t371 + t250 * t369) * t341; t251 * t386 + t252 * t387 + t253 * t388 + (-t254 * t382 - t255 * t381 - t256 * t380) * t392, t251 * t389 + t252 * t390 + t253 * t391 + (t254 * t379 + t255 * t378 + t256 * t377) * t392, t251 * t361 + t252 * t359 + t253 * t357 + m(4) + (t254 * t373 + t255 * t371 + t256 * t369) * t341;];
MX  = t1;
