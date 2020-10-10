% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:26
% EndTime: 2020-08-06 18:21:27
% DurationCPUTime: 1.33s
% Computational Cost: add. (4494->199), mult. (5745->333), div. (1308->11), fcn. (5232->18), ass. (0->145)
t529 = (m(3) * pkin(5));
t528 = 2 * qJ(2,1);
t527 = 2 * qJ(2,2);
t526 = 2 * qJ(2,3);
t448 = sin(qJ(3,3));
t431 = 0.1e1 / t448;
t454 = cos(qJ(3,3));
t513 = t431 * t454;
t432 = 0.1e1 / t448 ^ 2;
t460 = xDP(3);
t507 = t460 ^ 2 / pkin(3) ^ 2;
t491 = t432 * t507;
t450 = sin(qJ(3,2));
t434 = 0.1e1 / t450;
t456 = cos(qJ(3,2));
t511 = t434 * t456;
t435 = 0.1e1 / t450 ^ 2;
t489 = t435 * t507;
t452 = sin(qJ(3,1));
t437 = 0.1e1 / t452;
t458 = cos(qJ(3,1));
t509 = t437 * t458;
t438 = 0.1e1 / t452 ^ 2;
t487 = t438 * t507;
t525 = -2 * mrSges(3,1);
t524 = -2 * mrSges(2,3);
t523 = -2 * Ifges(3,4);
t447 = (Ifges(3,1) - Ifges(3,2));
t522 = 2 * t447;
t463 = (pkin(1) + pkin(5));
t520 = (mrSges(3,3) - mrSges(2,2));
t521 = (t520 + t529) * pkin(1);
t444 = legFrame(3,3);
t422 = sin(t444);
t425 = cos(t444);
t449 = sin(qJ(1,3));
t455 = cos(qJ(1,3));
t399 = -t422 * t449 + t425 * t455;
t402 = t422 * t455 + t425 * t449;
t412 = pkin(3) * t448 + qJ(2,3);
t409 = 0.1e1 / t412;
t461 = xDP(2);
t462 = xDP(1);
t381 = (t399 * t462 + t402 * t461) * t409;
t519 = t381 * t409;
t429 = pkin(6) + t463;
t518 = t381 * t429;
t445 = legFrame(2,3);
t423 = sin(t445);
t426 = cos(t445);
t451 = sin(qJ(1,2));
t457 = cos(qJ(1,2));
t400 = -t423 * t451 + t426 * t457;
t403 = t423 * t457 + t426 * t451;
t413 = pkin(3) * t450 + qJ(2,2);
t410 = 0.1e1 / t413;
t382 = (t400 * t462 + t403 * t461) * t410;
t517 = t382 * t410;
t516 = t382 * t429;
t446 = legFrame(1,3);
t424 = sin(t446);
t427 = cos(t446);
t453 = sin(qJ(1,1));
t459 = cos(qJ(1,1));
t401 = -t424 * t453 + t427 * t459;
t404 = t424 * t459 + t427 * t453;
t414 = pkin(3) * t452 + qJ(2,1);
t411 = 0.1e1 / t414;
t383 = (t401 * t462 + t404 * t461) * t411;
t515 = t383 * t411;
t514 = t383 * t429;
t512 = t431 * t460;
t510 = t434 * t460;
t508 = t437 * t460;
t506 = t448 * qJ(2,3);
t505 = t450 * qJ(2,2);
t504 = t452 * qJ(2,1);
t470 = 0.1e1 / pkin(3);
t503 = t460 * t470;
t492 = t454 * t512;
t393 = t412 * t449 + t429 * t455;
t396 = -t412 * t455 + t429 * t449;
t384 = t393 * t425 - t396 * t422;
t387 = t393 * t422 + t396 * t425;
t496 = (t384 * t462 + t387 * t461) * t409;
t366 = t492 + t496;
t440 = t454 ^ 2;
t465 = qJ(2,3) ^ 2;
t469 = pkin(3) ^ 2;
t480 = pkin(3) * t429 * t503;
t493 = -t429 ^ 2 - t469;
t357 = ((t366 * t429 - t480 * t513) * t448 + ((t440 * t469 - t465 + t493) * t448 + (t440 - 0.1e1) * pkin(3) * t526) * t381) * t431 * t519 + (t366 * t518 - ((t454 * t518 + t512) * t448 + qJ(2,3) * t431 * t503) * t432 * t460) * t409;
t360 = (t366 - t492 + t496 - t518) * t519;
t418 = mrSges(3,2) * t463 - Ifges(3,6);
t419 = mrSges(3,1) * t463 - Ifges(3,5);
t390 = -t418 * t448 + t419 * t454;
t405 = -mrSges(3,1) * t454 + mrSges(3,2) * t448;
t408 = m(2) * pkin(1) + m(3) * t463 + t520;
t464 = m(2) + m(3);
t415 = qJ(2,3) * t464 + mrSges(2,3);
t479 = -t464 * pkin(1) ^ 2 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + (-2 * mrSges(3,3) - t529) * pkin(5);
t483 = -mrSges(3,1) * t448 - mrSges(3,2) * t454;
t486 = t491 * t513;
t502 = ((qJ(2,3) * t524) - t447 * t440 - (t464 * t465) + t506 * t525 + t479) * t360 + t408 * t357 + t390 * t486 + 0.2e1 * (-(mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t448) * t454 - t521) * t360 + (t418 * t454 + t419 * t448) * t491 + (0.2e1 * (t415 - t483) * t366 + (t454 * t522 + (t405 * t526 + (0.4e1 * t440 - 0.2e1) * Ifges(3,4)) * t431) * t503) * t381;
t490 = t456 * t510;
t394 = t413 * t451 + t429 * t457;
t397 = -t413 * t457 + t429 * t451;
t385 = t394 * t426 - t397 * t423;
t388 = t394 * t423 + t397 * t426;
t495 = (t385 * t462 + t388 * t461) * t410;
t367 = t490 + t495;
t441 = t456 ^ 2;
t466 = qJ(2,2) ^ 2;
t358 = ((t367 * t429 - t480 * t511) * t450 + ((t441 * t469 - t466 + t493) * t450 + (t441 - 0.1e1) * pkin(3) * t527) * t382) * t434 * t517 + (t367 * t516 - ((t456 * t516 + t510) * t450 + qJ(2,2) * t434 * t503) * t435 * t460) * t410;
t361 = (t367 - t490 + t495 - t516) * t517;
t391 = -t418 * t450 + t419 * t456;
t406 = -mrSges(3,1) * t456 + mrSges(3,2) * t450;
t416 = qJ(2,2) * t464 + mrSges(2,3);
t482 = -mrSges(3,1) * t450 - mrSges(3,2) * t456;
t485 = t489 * t511;
t501 = ((qJ(2,2) * t524) - t447 * t441 - (t464 * t466) + t505 * t525 + t479) * t361 + t408 * t358 + t391 * t485 + 0.2e1 * (-(mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t450) * t456 - t521) * t361 + (t418 * t456 + t419 * t450) * t489 + (0.2e1 * (t416 - t482) * t367 + (t456 * t522 + (t406 * t527 + (0.4e1 * t441 - 0.2e1) * Ifges(3,4)) * t434) * t503) * t382;
t488 = t458 * t508;
t395 = t414 * t453 + t429 * t459;
t398 = -t414 * t459 + t429 * t453;
t386 = t395 * t427 - t398 * t424;
t389 = t395 * t424 + t398 * t427;
t494 = (t386 * t462 + t389 * t461) * t411;
t368 = t488 + t494;
t442 = t458 ^ 2;
t467 = qJ(2,1) ^ 2;
t359 = ((t368 * t429 - t480 * t509) * t452 + ((t442 * t469 - t467 + t493) * t452 + (t442 - 0.1e1) * pkin(3) * t528) * t383) * t437 * t515 + (t368 * t514 - ((t458 * t514 + t508) * t452 + qJ(2,1) * t437 * t503) * t438 * t460) * t411;
t362 = (t368 - t488 + t494 - t514) * t515;
t392 = -t418 * t452 + t419 * t458;
t407 = -mrSges(3,1) * t458 + mrSges(3,2) * t452;
t417 = qJ(2,1) * t464 + mrSges(2,3);
t481 = -mrSges(3,1) * t452 - mrSges(3,2) * t458;
t484 = t487 * t509;
t500 = ((qJ(2,1) * t524) - t447 * t442 - (t464 * t467) + t504 * t525 + t479) * t362 + t408 * t359 + t392 * t484 + 0.2e1 * (-(mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t452) * t458 - t521) * t362 + (t418 * t458 + t419 * t452) * t487 + (0.2e1 * (t417 - t481) * t368 + (t458 * t522 + (t407 * t528 + (0.4e1 * t442 - 0.2e1) * Ifges(3,4)) * t437) * t503) * t383;
t378 = t381 ^ 2;
t499 = -t357 * t464 + t360 * t408 + t405 * t486 - t378 * t415 + t483 * (t378 + t491);
t379 = t382 ^ 2;
t498 = -t358 * t464 + t361 * t408 + t406 * t485 - t379 * t416 + t482 * (t379 + t489);
t380 = t383 ^ 2;
t497 = -t359 * t464 + t362 * t408 + t407 * t484 - t380 * t417 + t481 * (t380 + t487);
t1 = [(t497 * t386 + t500 * t401) * t411 + (t498 * t385 + t501 * t400) * t410 + (t499 * t384 + t502 * t399) * t409; (t497 * t389 + t500 * t404) * t411 + (t498 * t388 + t501 * t403) * t410 + (t499 * t387 + t502 * t402) * t409; t497 * t509 + t498 * t511 + t499 * t513 + ((Ifges(3,3) * t484 - t359 * t407 - t362 * t392 + (t442 * t523 + ((mrSges(3,1) * qJ(2,1)) - t447 * t452) * t458 - mrSges(3,2) * t504 + Ifges(3,4)) * t380) * t437 + (Ifges(3,3) * t485 - t358 * t406 - t361 * t391 + (t441 * t523 + ((mrSges(3,1) * qJ(2,2)) - t447 * t450) * t456 - mrSges(3,2) * t505 + Ifges(3,4)) * t379) * t434 + (Ifges(3,3) * t486 - t357 * t405 - t360 * t390 + (t440 * t523 + ((mrSges(3,1) * qJ(2,3)) - t447 * t448) * t454 - mrSges(3,2) * t506 + Ifges(3,4)) * t378) * t431) * t470;];
taucX  = t1;
