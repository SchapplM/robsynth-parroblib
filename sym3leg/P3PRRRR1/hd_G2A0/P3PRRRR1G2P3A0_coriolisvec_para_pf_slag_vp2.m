% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G2P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:28
% EndTime: 2020-03-09 21:16:29
% DurationCPUTime: 0.92s
% Computational Cost: add. (1023->124), mult. (2925->262), div. (2598->10), fcn. (3726->18), ass. (0->121)
t453 = sin(qJ(2,1));
t433 = 0.1e1 / t453;
t458 = cos(qJ(3,1));
t441 = 0.1e1 / t458;
t460 = xDP(3);
t463 = 0.1e1 / pkin(2);
t447 = legFrame(1,2);
t426 = sin(t447);
t429 = cos(t447);
t461 = xDP(2);
t462 = xDP(1);
t464 = t426 * t461 - t429 * t462;
t442 = 0.1e1 / t458 ^ 2;
t452 = sin(qJ(3,1));
t459 = cos(qJ(2,1));
t485 = t442 * t452 * t459;
t399 = (-t441 * t460 - t464 * t485) * t463 * t433;
t405 = t464 * t463 * t441;
t491 = t453 * t458;
t500 = t405 * t452;
t506 = t399 * t459;
t507 = t399 * t452;
t384 = ((-t453 * t500 + t458 * t506) * t441 * t399 + (t459 * t405 - t491 * t507) * t442 * t405) * t433;
t396 = t399 ^ 2;
t402 = t405 ^ 2;
t503 = t402 * t441;
t390 = (-t396 * t458 - t503) * t433 * pkin(2);
t420 = -Ifges(3,5) * t452 - Ifges(3,6) * t458;
t443 = Ifges(3,1) - Ifges(3,2);
t444 = mrSges(2,2) - mrSges(3,3);
t470 = -t458 * mrSges(3,1) + t452 * mrSges(3,2);
t408 = -(mrSges(2,1) - t470) * t459 + t453 * t444;
t440 = t458 ^ 2;
t488 = t452 * t503;
t492 = t452 * t458;
t512 = -Ifges(3,1) - Ifges(2,3);
t514 = -Ifges(3,6) / 0.2e1;
t515 = Ifges(3,5) / 0.2e1;
t516 = 0.2e1 * t440;
t519 = -2 * Ifges(3,4);
t467 = t433 * (0.2e1 * ((t405 * t515 + t443 * t507) * t458 + t500 * t514 + (t516 - 0.1e1) * t399 * Ifges(3,4)) * t405 + t408 * t390 + (t443 * t440 + t492 * t519 + t512) * t384 - t420 * t488);
t423 = mrSges(3,1) * t452 + mrSges(3,2) * t458;
t497 = t423 * t453;
t525 = t467 * t485 + (-Ifges(3,3) * t488 - t420 * t384 - t390 * t497 + t396 * (Ifges(3,4) * t516 + t443 * t492 - Ifges(3,4))) * t441;
t451 = sin(qJ(2,2));
t432 = 0.1e1 / t451;
t456 = cos(qJ(3,2));
t438 = 0.1e1 / t456;
t446 = legFrame(2,2);
t425 = sin(t446);
t428 = cos(t446);
t465 = t425 * t461 - t428 * t462;
t439 = 0.1e1 / t456 ^ 2;
t450 = sin(qJ(3,2));
t457 = cos(qJ(2,2));
t486 = t439 * t450 * t457;
t398 = (-t438 * t460 - t465 * t486) * t463 * t432;
t404 = t465 * t463 * t438;
t493 = t451 * t456;
t501 = t404 * t450;
t508 = t398 * t457;
t509 = t398 * t450;
t383 = ((-t451 * t501 + t456 * t508) * t438 * t398 + (t457 * t404 - t493 * t509) * t439 * t404) * t432;
t395 = t398 ^ 2;
t401 = t404 ^ 2;
t504 = t401 * t438;
t389 = (-t395 * t456 - t504) * t432 * pkin(2);
t419 = -Ifges(3,5) * t450 - Ifges(3,6) * t456;
t471 = -t456 * mrSges(3,1) + t450 * mrSges(3,2);
t407 = -(mrSges(2,1) - t471) * t457 + t451 * t444;
t437 = t456 ^ 2;
t489 = t450 * t504;
t494 = t450 * t456;
t517 = 0.2e1 * t437;
t468 = t432 * (0.2e1 * ((t404 * t515 + t443 * t509) * t456 + t501 * t514 + (t517 - 0.1e1) * t398 * Ifges(3,4)) * t404 + t407 * t389 + (t443 * t437 + t494 * t519 + t512) * t383 - t419 * t489);
t422 = mrSges(3,1) * t450 + mrSges(3,2) * t456;
t498 = t422 * t451;
t524 = t468 * t486 + (-Ifges(3,3) * t489 - t419 * t383 - t389 * t498 + t395 * (Ifges(3,4) * t517 + t443 * t494 - Ifges(3,4))) * t438;
t449 = sin(qJ(2,3));
t431 = 0.1e1 / t449;
t454 = cos(qJ(3,3));
t435 = 0.1e1 / t454;
t445 = legFrame(3,2);
t424 = sin(t445);
t427 = cos(t445);
t466 = t424 * t461 - t427 * t462;
t436 = 0.1e1 / t454 ^ 2;
t448 = sin(qJ(3,3));
t455 = cos(qJ(2,3));
t487 = t436 * t448 * t455;
t397 = (-t435 * t460 - t466 * t487) * t463 * t431;
t403 = t466 * t463 * t435;
t495 = t449 * t454;
t502 = t403 * t448;
t510 = t397 * t455;
t511 = t397 * t448;
t382 = ((-t449 * t502 + t454 * t510) * t435 * t397 + (t455 * t403 - t495 * t511) * t436 * t403) * t431;
t394 = t397 ^ 2;
t400 = t403 ^ 2;
t505 = t400 * t435;
t388 = (-t394 * t454 - t505) * t431 * pkin(2);
t418 = -Ifges(3,5) * t448 - Ifges(3,6) * t454;
t472 = -t454 * mrSges(3,1) + t448 * mrSges(3,2);
t406 = -(mrSges(2,1) - t472) * t455 + t449 * t444;
t434 = t454 ^ 2;
t490 = t448 * t505;
t496 = t448 * t454;
t518 = 0.2e1 * t434;
t469 = t431 * (0.2e1 * ((t403 * t515 + t443 * t511) * t454 + t502 * t514 + (t518 - 0.1e1) * t397 * Ifges(3,4)) * t403 + t406 * t388 + (t443 * t434 + t496 * t519 + t512) * t382 - t418 * t490);
t421 = mrSges(3,1) * t448 + mrSges(3,2) * t454;
t499 = t421 * t449;
t523 = t469 * t487 + (-Ifges(3,3) * t490 - t418 * t382 - t388 * t499 + t394 * (Ifges(3,4) * t518 + t443 * t496 - Ifges(3,4))) * t435;
t513 = t444 / 0.2e1;
t430 = -m(1) - m(2) - m(3);
t484 = t431 * (t406 * t382 + t430 * t388 - t490 * t499 + (-mrSges(2,1) * t394 + t472 * (t394 + t400)) * t449 - 0.2e1 * (t397 * t513 + t421 * t403) * t510);
t483 = t432 * (t407 * t383 + t430 * t389 - t489 * t498 + (-mrSges(2,1) * t395 + t471 * (t395 + t401)) * t451 - 0.2e1 * (t398 * t513 + t422 * t404) * t508);
t482 = t433 * (t408 * t384 + t430 * t390 - t488 * t497 + (-mrSges(2,1) * t396 + t470 * (t396 + t402)) * t453 - 0.2e1 * (t399 * t513 + t423 * t405) * t506);
t475 = t435 * t484;
t474 = t438 * t483;
t473 = t441 * t482;
t1 = [(t426 * t491 - t429 * t452) * t473 + (t425 * t493 - t428 * t450) * t474 + (t424 * t495 - t427 * t448) * t475 + (t523 * t427 + t524 * t428 + t525 * t429) * t463; (t426 * t452 + t429 * t491) * t473 + (t425 * t450 + t428 * t493) * t474 + (t424 * t448 + t427 * t495) * t475 + (-t523 * t424 - t524 * t425 - t525 * t426) * t463; t459 * t482 + t457 * t483 + t455 * t484 + (-t435 * t469 - t438 * t468 - t441 * t467) * t463;];
taucX  = t1;
