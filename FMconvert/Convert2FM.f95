

!  This program tries to automate most of the conversion of a program to use the FM package
!  for multiple precision computation.

!  To convert everything would require almost as much syntax analysis as a full compiler.
!  Instead, the most commonly used Fortran features are handled, by using some shortcuts
!  and simplifying assumptions about the syntax.

!  C2FM.INP is the input file  -- the non multiple precision program.
!  C2FM.OUT is the output file -- the FM version.
!  Several other files are used for temporary information, these can be ignored.

!  Because not all the variables in the original program may need to be multiple precision,
!  this conversion program uses the convention that the constant in the KIND part of the
!  type declaration statement is a code telling when to convert a variable to multiple precision.
!  Make a copy of the original program and call it C2FM.INP, then change each declaration there
!  for which the variables will be multiple precision in the new version to one of the following
!  three types.  Use exactly these types, with all upper case letters and one blank.

!  REAL (KIND(3.1D1))     will be converted to  TYPE (FM),
!  COMPLEX (KIND(3.1D1))  will be converted to  TYPE (ZM),
!  INTEGER (KIND(31))     will be converted to  TYPE (IM).

!  For an existing program, decide which variables need to become multiple precision, and
!  change each declaration to one of these.  In the common case of a program that will not
!  need multiple precision complex or integers, and all reals are double precision before
!  the conversion, this can usually be done with one global change in an editor.

!  Calls to functions DBLE, FLOAT, REAL, SNGL that are in lines that reference multiple precision
!  variables or routines will be converted to calls to TO_FM.  This could give error messages if
!  a REAL call in the original program uses the optional KIND argument, and those will have to
!  be fixed by hand.

!  Any multiple precision functions defined in your program should be declared EXTERNAL along
!  with the proper kind code.  Omitting the external attribute might result in a compiler error
!  message because of trying to initialize it.

!  If your program uses variables that will become multiple precision in a complicated way, then
!  Convert2FM is likely to miss some changes.  These will have to be found and changed by hand,
!  but most of the time they will cause compiler error messages that show where the missed
!  changes were.  For example, if your program has multiple precision components of derived types
!  or pointer aliases to multiple precision variables, Convert2FM might miss some changes.

!  The most time-consuming parts of conversion are usually changing constants,
!     X = Y/3.7   should become   X = Y/TO_FM('3.7'),
!  changing type declarations,
!     REAL (KIND(1.0D0)) :: X, Y   should become   TYPE (FM), SAVE :: X, Y
!  and inserting FM interface instructions like
!     USE FMZM
!     CALL FM_ENTER_USER_ROUTINE
!     CALL FM_EXIT_USER_ROUTINE.
!  Convert2FM should be able to do these automatically.

!  Converting read/write statements is done by converting between machine precision and multiple
!  precision, so that any associated formats do not have to be changed.  These will often need
!  to be changed by hand to get the desired multiple precision formats.  This is especially true
!  for TYPE (IM) numbers, since the machine precision integer overflow threshold is very low,
!  so writing TO_INT(K) is seldom correct.

!  For example,
!      WRITE (KW,120) N,H,TOL
!  becomes
!      WRITE (KW,120) N,TO_DP(H),TO_DP(TOL)
!  after Convert2FM, which writes the multiple precision values H and TOL only to double
!  precision accuracy.

!  Reads are similar but messier:
!      READ (KR,*) (X(L),L=1,N)
!  becomes
!      ALLOCATE( C2FM_TEMP(1:N) )
!      READ (KR,*) (C2FM_TEMP(L),L=1,N)
!      DO L = 1, N
!         X(L) = C2FM_TEMP(L)
!      ENDDO
!      DEALLOCATE(C2FM_TEMP)
!  where X is type (fm) in the new version.  C2FM_TEMP is a double precision variable in module
!  C2FM_READS that Convert2FM puts at the start of the new program.


      MODULE C2FM_VARS


!  These are the global variables for C2FM.

!  LINE contains one subprogram at a time while the input file is analyzed and converted to FM.

!  ROUTINE_NAME     is the name of the current routine.

!  ROUTINE_TYPE = 1 for the main program
!               = 2 for subroutines
!               = 3 for functions
!               = 4 for modules

!  N_ARGS is the number of arguments in the current routine.

!  ROUTINE_ARGUMENT(K) is the variable name of the Kth argument in the current routine.

!  ARGUMENT_TYPE(K) = 0 if that variable will not become multiple precision
!                   = 1 for TYPE (FM)
!                   = 2 for TYPE (ZM)
!                   = 3 for TYPE (IM)

!  N_VARS is the number of variables in the current routine.

!  ROUTINE_VARIABLE(K) is the variable name of the Kth variable in the current routine that
!                      will become multiple precision.

!  VARIABLE_INITIALIZATION(K) gives any initialization expression in the declaration of
!                             the original program.  These must be made into executable
!                             statements and put at the top of the routine.

!  VARIABLE_LOCATION(K) records whether that variable is an input argument, a module variable,
!                       or a local variable.

!  ARRAY_SIZE(K) gives the size of a multiple precision variable array.

!  VARIABLE_TYPE(K) = 0, 1, 2, 3 gives these local variable types as with ARGUMENT_TYPE.

!  The next five arrays form a database showing which routine input arguments are multiple
!  precision types.

!  N_ROUTINES is the number of routines in the program.

!  R_NAME(J) is the routine name of the Jth subprogram in the sorted list of names.

!  R_TYPE(J) is the routine type of the Jth subprogram.

!  R_N_ARGS(J) is the number of arguments for the Jth routine.

!  R_ARG_START(J) gives the starting location of the list of argument types for routine J.

!  R_ARGS( R_ARG_START(J) : R_ARG_START(J) + R_N_ARGS(J) - 1 ) gives the argument type of each
!          input argument to the Jth routine, using the same codes as ARGUMENT_TYPE.

!  The next five arrays form a database showing which modules are used.

!  R_N_USED(J) is the number of USE statements in the Jth routine.

!  R_USED_START(J) gives the starting location of the list of used modules for routine J.

!  R_USED( R_USED_START(J) : R_USED_START(J) + R_N_USED(J) - 1 ) gives the module routine names
!          used in the Jth routine.

!  R_N_ALL_USED(J) is the number of modules used in the Jth routine, including indirect usage
!          where a directly used module uses another module, etc.

!  R_ALL_USED( J , 1 : R_N_ALL_USED(J) ) gives the module routine numbers directly and indirectly
!          used in the Jth routine.

!  AFTER_READ holds new conversion statements used after the new read statement gets double
!          precision inputs and then these statements convert those values to multiple precision.


      CHARACTER(132),  SAVE :: LINE(10000)
      CHARACTER(32),   SAVE :: ROUTINE_NAME,ROUTINE_ARGUMENT(100),ROUTINE_VARIABLE(1000)
      CHARACTER(32),   SAVE :: R_NAME(1000),R_USED(10000)
      CHARACTER(132),  SAVE :: VARIABLE_INITIALIZATION(1000)
      CHARACTER(9999), SAVE :: TEMP_STATEMENT,TEMP_ST
      CHARACTER(8),    SAVE :: VARIABLE_LOCATION(1000)
      CHARACTER(132),  SAVE :: AFTER_READ(100)
      INTEGER, SAVE :: ARGUMENT_TYPE(100),VARIABLE_TYPE(1000),ARRAY_SIZE(1000),ROUTINE_TYPE,     &
                       N_ARGS,N_VARS,QUOTE_LEVEL(9999)
      INTEGER, SAVE :: N_ROUTINES,R_TYPE(1000),R_N_ARGS(1000),R_ARG_START(1000),R_ARGS(100000),  &
                       R_N_USED(1000),R_USED_START(1000),R_N_ALL_USED(1000),N_AFTER_READ,        &
                       N_AFTER_ALLOC
      INTEGER, SAVE, DIMENSION(:,:), ALLOCATABLE :: R_ALL_USED

      END MODULE C2FM_VARS




      PROGRAM C2FM

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: FIRST,J,LAST
      LOGICAL :: FILE_ENDED

      OPEN(22,FILE='C2FM.INP')
      OPEN(23,FILE='C2FM.OUT')
      OPEN(31,FILE='C2FM.ARGS')
      OPEN(32,FILE='C2FM.VARS')
      OPEN(33,FILE='C2FM.MODS')
      OPEN(34,FILE='C2FM.PASS1')
      OPEN(35,FILE='C2FM.VARS2')
      OPEN(36,FILE='C2FM.PASS2')

!             Not all programs use modules.  Write a blank MODS file to erase any information
!             left over from previous conversions.

      WRITE (33,*) ' '
      CLOSE(33)
      OPEN(33,FILE='C2FM.MODS')

!             Make a first pass over the program to find out which arguments to each subprogram
!             will be multiple precision.

      CALL GET_ARGUMENTS

!             Sort the R_NAME database by routine names to speed up later searches.

      CALL SORT_NAMES

!             Expand the R_NAME database by adding any modules used by other modules to
!             the list of used modules in each routine.

      CALL USE_CHECK

!             Make a second variable list for each routine, including any module variables
!             that are used in each routine.

      CALL VAR2_LIST

!             The second pass adds code to initialize the multiple precision variables and
!             record entry and exit to routines so that multiple precision temporary variables
!             do not get discarded too soon.

      REWIND(31)
      REWIND(32)
      REWIND(33)
      REWIND(34)
      REWIND(35)

      DO J = 1, 10000

!             Read the next routine.

         CALL READ_ROUTINE(34,FIRST,LAST,FILE_ENDED)
         IF (FILE_ENDED .AND. LAST <= 0) EXIT

!             Get the information about local multiple precision variables.

         CALL LOCAL_VARS(32)

!             Write the next routine to the output file.

         CALL WRITE_PASS2(LAST)

         IF (FILE_ENDED) EXIT

      ENDDO

!             The third pass converts floating point constants that appear in expressions
!             involving multiple precision variables, and also converts read, write, print,
!             and allocate statements.

      REWIND(31)
      REWIND(32)
      REWIND(33)
      REWIND(34)
      REWIND(35)
      REWIND(36)

!             Put a module with two allocatable arrays at the top of the FM version of the program.
!             These are used by the new code that handles read statements.

      WRITE (23,"(A)") ' '
      WRITE (23,"(A)") '      MODULE C2FM_READS'
      WRITE (23,"(A)") '         REAL (KIND(1.0D0)), DIMENSION(:),   ALLOCATABLE :: C2FM_TEMP'
      WRITE (23,"(A)") '         REAL (KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: C2FM_TEMP2'
      WRITE (23,"(A)") '      END MODULE C2FM_READS'
      WRITE (23,"(A)") ' '

      DO J = 1, 10000

!             Read the next routine.

         CALL READ_ROUTINE(36,FIRST,LAST,FILE_ENDED)
         IF (FILE_ENDED .AND. LAST <= 0) STOP

!             Get the information about local multiple precision variables.

         CALL LOCAL_VARS(35)

!             Write the next routine to the output file.

         CALL WRITE_PASS3(LAST)

         IF (FILE_ENDED) STOP

      ENDDO

      END PROGRAM C2FM

      SUBROUTINE ALL_CAPS(STRING,N,RESULT)

!  Convert STRING(1:N) to all upper case and return it an RESULT(1:N).

      IMPLICIT NONE
      INTEGER :: N
      CHARACTER(N) :: STRING,RESULT
      CHARACTER(26), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz',  &
                                  UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      INTEGER :: K,KA

      RESULT(1:N) = STRING(1:N)
      DO K = 1, N
         KA = INDEX(LOWER_CASE,STRING(K:K))
         IF (KA > 0) RESULT(K:K) = UPPER_CASE(KA:KA)
      ENDDO

      END SUBROUTINE ALL_CAPS

      SUBROUTINE CONVERT_CONSTANTS(START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT,J_START,J_END)

!  TEMP_STATEMENT contains the statement.
!  Look for an executable statement where constants should be converted to
!  multiple precision, and insert calls to TO_FM, etc.
!
!  START_OF_STATEMENT is the location of the first nonblank character in the statement.
!  LAST_NONBLANK is the location of the last nonblank character in the statement,
!                including any trailing comments after !.
!  END_OF_STATEMENT is the location of the last nonblank character in the non-comment
!                   part of the statement.
!  J_START and J_END sometimes limit the columns where constants will be changed.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(63), PARAMETER :: NAME_CHARS =  &
        'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789'
      CHARACTER(32), SAVE :: VAR_NAME,NAME_OF_ROUTINE
      INTEGER :: END_OF_STATEMENT,EXP_CHAR,START_OF_STATEMENT,J,J1,JA,JB,J_START,J_END,K,KA,KB,  &
                 K_END,L,L2,NARG,LAST_NONBLANK,SKIP_UNTIL,LEVEL_P,LEVEL_Q1,LEVEL_Q2,INDEX_ULC
      LOGICAL :: IN_NAME,IN_STRING,IN_CONSTANT,IN_CALL,INTEGER_CONSTANT,  &
                 DIGIT_IN_CONSTANT,ROUTINE_CALL,FM_VAR_FOUND,PREVIOUS_CHAR_NAME

!             See if any multiple precision variable names appear in this statement.
!             Also look for calls to functions or subroutines.

      ROUTINE_CALL = .FALSE.
      FM_VAR_FOUND = .FALSE.
      IN_NAME = .FALSE.
      K_END = J_END
      DO J = START_OF_STATEMENT, END_OF_STATEMENT+1
         IF (.NOT. IN_NAME) THEN
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 > 0 .AND. J1 < 53) THEN
                 IN_NAME = .TRUE.
                 JA = J
             ENDIF
         ELSE
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 == 0) THEN
                 IN_NAME = .FALSE.
                 JB = J - 1
                 VAR_NAME = ' '
                 CALL ALL_CAPS(TEMP_STATEMENT(JA:JB),JB-JA+1,VAR_NAME)
                 DO K = 1, N_VARS
                    IF (VAR_NAME == ROUTINE_VARIABLE(K)) THEN
                        FM_VAR_FOUND = .TRUE.
                        EXIT
                    ENDIF
                 ENDDO
                 DO K = 1, N_ROUTINES
                    IF (VAR_NAME == R_NAME(K)) THEN
                        ROUTINE_CALL = .TRUE.
                        NAME_OF_ROUTINE = VAR_NAME
                        EXIT
                    ENDIF
                 ENDDO
             ENDIF
         ENDIF
      ENDDO

      IF (.NOT.(FM_VAR_FOUND .OR. ROUTINE_CALL)) RETURN

!             Look for any calls to DBLE, FLOAT, REAL, SNGL.
!             Convert them to TO_FM.

      IN_NAME = .FALSE.
      J = START_OF_STATEMENT
  110    IF (.NOT. IN_NAME) THEN
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 > 0 .AND. J1 < 53) THEN
                 IN_NAME = .TRUE.
                 JA = J
             ENDIF
         ELSE
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 == 0) THEN
                 IN_NAME = .FALSE.
                 JB = J - 1
                 VAR_NAME = ' '
                 CALL ALL_CAPS(TEMP_STATEMENT(JA:JB),JB-JA+1,VAR_NAME)
                 IF (VAR_NAME == 'DBLE' .OR. VAR_NAME == 'REAL' .OR. VAR_NAME == 'SNGL') THEN
                     DO K = LAST_NONBLANK, JB+1, -1
                        TEMP_STATEMENT(K+1:K+1) = TEMP_STATEMENT(K:K)
                     ENDDO
                     TEMP_STATEMENT(JA:JB+1) = 'TO_FM'
                     LAST_NONBLANK = LAST_NONBLANK + 1
                     END_OF_STATEMENT = END_OF_STATEMENT + 1
                 ENDIF
                 IF (VAR_NAME == 'FLOAT') THEN
                     TEMP_STATEMENT(JA:JB) = 'TO_FM'
                 ENDIF
             ENDIF
         ENDIF
      J = J + 1
      IF (J <= END_OF_STATEMENT+1) GO TO 110

!             Scan for constants to convert.
!             All non-integer constants are converted, except in format statements.

      IN_NAME = .FALSE.
      IN_STRING = .FALSE.
      IN_CONSTANT = .FALSE.
      IN_CALL = .FALSE.
      SKIP_UNTIL = 0
      J = START_OF_STATEMENT
  120    IF (SKIP_UNTIL > 0 .AND. J < SKIP_UNTIL) GO TO 140
         IF (J < J_START .OR. J > K_END) GO TO 140
         IF (.NOT. (IN_NAME .OR. IN_STRING .OR. IN_CONSTANT)) THEN
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 > 0 .AND. J1 < 53) THEN
                 IN_NAME = .TRUE.
                 JA = J
             ELSE IF (J1 > 53 .OR. TEMP_STATEMENT(J:J) == '.') THEN
                 IN_CONSTANT = .TRUE.
                 JA = J
                 EXP_CHAR = 0
                 IF (J1 > 53) THEN
                     INTEGER_CONSTANT = .TRUE.
                     DIGIT_IN_CONSTANT = .TRUE.
                 ELSE
                     INTEGER_CONSTANT = .FALSE.
                     DIGIT_IN_CONSTANT = .FALSE.
                 ENDIF
             ELSE IF (TEMP_STATEMENT(J:J) == "'") THEN
                 IN_STRING = .TRUE.
                 JA = J
             ELSE IF (TEMP_STATEMENT(J:J) == '"') THEN
                 IN_STRING = .TRUE.
                 JA = J
             ENDIF
         ELSE IF (IN_NAME) THEN
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 == 0) THEN
                 IN_NAME = .FALSE.
                 JB = J - 1
                 VAR_NAME = TEMP_STATEMENT(JA:JB)
                 IF (INDEX_ULC(VAR_NAME,'FORMAT') == 1 .OR. INDEX_ULC(VAR_NAME,'READ')  == 1 .OR.  &
                     INDEX_ULC(VAR_NAME,'WRITE')  == 1 .OR. INDEX_ULC(VAR_NAME,'PRINT') == 1) THEN
                     GO TO 150
                 ENDIF
                 DO K = 1, N_ROUTINES
                    IF (VAR_NAME == R_NAME(K)) THEN
                        ROUTINE_CALL = .TRUE.
                        NAME_OF_ROUTINE = VAR_NAME
                        IN_CALL = .TRUE.
                        EXIT
                    ENDIF
                 ENDDO
             ENDIF
         ELSE IF (IN_STRING) THEN
             IF (TEMP_STATEMENT(J:J) == TEMP_STATEMENT(JA:JA)) THEN
                 IN_STRING = .FALSE.
             ENDIF
         ELSE IF (IN_CONSTANT) THEN
             J1 = INDEX(NAME_CHARS,TEMP_STATEMENT(J:J))
             IF (J1 > 53) THEN
                 DIGIT_IN_CONSTANT = .TRUE.
                 IF (J >= END_OF_STATEMENT) THEN
                     J = J + 1
                     GO TO 130
                 ENDIF
                 GO TO 140
             ENDIF
             IF (TEMP_STATEMENT(J:J) == '.') THEN
                 INTEGER_CONSTANT = .FALSE.
                 IF (J >= END_OF_STATEMENT) THEN
                     J = J + 1
                     GO TO 130
                 ENDIF
                 GO TO 140
             ENDIF
             IF (EXP_CHAR == 0) THEN
                 IF (TEMP_STATEMENT(J:J) == 'E' .OR. TEMP_STATEMENT(J:J) == 'e' .OR.  &
                     TEMP_STATEMENT(J:J) == 'D' .OR. TEMP_STATEMENT(J:J) == 'd') THEN
                     EXP_CHAR = J
                     INTEGER_CONSTANT = .FALSE.
                     IF (J >= END_OF_STATEMENT) THEN
                         J = J + 1
                         GO TO 130
                     ENDIF
                     GO TO 140
                 ENDIF
             ENDIF
             IF ((TEMP_STATEMENT(J:J) == '+' .OR. TEMP_STATEMENT(J:J) == '-') .AND.  &
                 EXP_CHAR > 0 .AND. J == EXP_CHAR+1) THEN
                 IF (J >= END_OF_STATEMENT) THEN
                     J = J + 1
                     GO TO 130
                 ENDIF
                 GO TO 140
             ENDIF

  130        JB = J - 1
             IN_CONSTANT = .FALSE.
             IF (.NOT. DIGIT_IN_CONSTANT) GO TO 140

!             End of constant found, expand the line and convert it.

             IF (.NOT. INTEGER_CONSTANT) THEN
                 DO K = LAST_NONBLANK, JB+1, -1
                    TEMP_STATEMENT(K+9:K+9) = TEMP_STATEMENT(K:K)
                 ENDDO
                 TEMP_STATEMENT(JB+8:JB+9) = "')"
                 DO K = JB, JA, -1
                    TEMP_STATEMENT(K+7:K+7) = TEMP_STATEMENT(K:K)
                 ENDDO
                 TEMP_STATEMENT(JA:JA+6) = "TO_FM('"
                 SKIP_UNTIL = JB + 10
                 K_END = MIN(9999,K_END+9)
                 LAST_NONBLANK = LAST_NONBLANK + 9
                 END_OF_STATEMENT = END_OF_STATEMENT + 9
             ENDIF
         ENDIF
  140 J = J + 1
      IF (J <= END_OF_STATEMENT+1) GO TO 120

!             Convert complex constants.

!             At this point, a complex constant like  (1.23,4.56)  will have been converted
!             to  (TO_FM('1.23'),TO_FM('4.56'))  if any multiple precision variables also
!             appear in this statement.  Now this needs to become
!             CMPLX(TO_FM('1.23'),TO_FM('4.56')) so it will be converted to type ZM.

  150 PREVIOUS_CHAR_NAME = .FALSE.
      SKIP_UNTIL = 0
      J1 = 0
      J = 1

  160    IF (SKIP_UNTIL > 0 .AND. J < SKIP_UNTIL) GO TO 170
         IF (J < J_START .OR. J > K_END) GO TO 170
         IF (TEMP_STATEMENT(J:J) == ' ') GO TO 170
         IF (.NOT. PREVIOUS_CHAR_NAME .AND. J1 == 0) THEN
             IF (TEMP_STATEMENT(J:J) == '(') THEN
                 J1 = 1
                 JA = J
                 GO TO 170
             ENDIF
         ENDIF
         IF (J1 == 1 .OR. J1 == 3) THEN
             IF (INDEX_ULC(TEMP_STATEMENT(J:J+5),'TO_FM(') == 1) THEN
                 J1 = J1 + 1
                 J = J + 6
                 GO TO 160
             ELSE
                 IF (TEMP_STATEMENT(J:J) == '+' .OR. TEMP_STATEMENT(J:J) == '-') THEN
                     GO TO 170
                 ELSE
                     J1 = 0
                     PREVIOUS_CHAR_NAME = .FALSE.
                     GO TO 160
                 ENDIF
             ENDIF
         ENDIF
         IF (J1 == 2) THEN
             IF (TEMP_STATEMENT(J:J) == ',') THEN
                 J1 = 3
                 GO TO 170
             ELSE
                 GO TO 170
             ENDIF
         ENDIF
         IF (J1 == 4 .OR. J1 == 5) THEN
             IF (TEMP_STATEMENT(J:J) == ')') THEN
                 J1 = J1 + 1
                 GO TO 170
             ELSE
                 GO TO 170
             ENDIF
         ENDIF

         IF (INDEX(NAME_CHARS,TEMP_STATEMENT(J:J)) > 0) THEN
             PREVIOUS_CHAR_NAME = .TRUE.
         ELSE
             PREVIOUS_CHAR_NAME = .FALSE.
         ENDIF

  170    IF (J1 == 6) THEN
             JB = J

!             End of complex constant found, expand the line and convert it.

             DO K = LAST_NONBLANK, JA, -1
                TEMP_STATEMENT(K+5:K+5) = TEMP_STATEMENT(K:K)
             ENDDO
             TEMP_STATEMENT(JA:JA+4) = "CMPLX"
             SKIP_UNTIL = JB + 6
             K_END = MIN(9999,K_END+5)
             LAST_NONBLANK = LAST_NONBLANK + 5
             END_OF_STATEMENT = END_OF_STATEMENT + 5

             J1 = 0
             PREVIOUS_CHAR_NAME = .FALSE.
         ENDIF

      J = J + 1
      IF (J <= END_OF_STATEMENT+1) GO TO 160

!             Now that constants in this statement have been converted to multiple precision,
!             look for any that are arguments to a function or subroutine and are not multiple
!             precision arguments.

      IF (INDEX_ULC(TEMP_STATEMENT(1:END_OF_STATEMENT),'TO_FM(') <= 0 .AND.  &
          INDEX_ULC(TEMP_STATEMENT(1:END_OF_STATEMENT),'TO_IM(') <= 0) RETURN

      DO J = 1, N_ROUTINES
         KA = INDEX_ULC(TEMP_STATEMENT(1:END_OF_STATEMENT),TRIM(R_NAME(J)))
  180    IF (KA <= 0) CYCLE
         KB = KA + LEN_TRIM(R_NAME(J)) - 1

!             Make sure the match is the full name and not a substring.

         IF (KA > 1) THEN
             IF (INDEX(NAME_CHARS,TEMP_STATEMENT(KA-1:KA-1)) > 0) CYCLE
         ENDIF
         IF (KB < END_OF_STATEMENT) THEN
             IF (INDEX(NAME_CHARS,TEMP_STATEMENT(KB+1:KB+1)) > 0) CYCLE
         ELSE
             RETURN
         ENDIF

!             Scan the argument list.

         LEVEL_P = 0
         LEVEL_Q1 = 0
         LEVEL_Q2 = 0
         NARG = 0
         JA = 0
         JB = -1
         DO K = KB+1, END_OF_STATEMENT
            IF (TEMP_STATEMENT(K:K) == '(' .AND. LEVEL_Q1 == 0 .AND. LEVEL_Q2 == 0) THEN
                LEVEL_P = LEVEL_P + 1
                IF (LEVEL_P == 1) THEN
                    JA = K + 1
                    NARG = 1
                ENDIF
            ENDIF
            IF (TEMP_STATEMENT(K:K) == ')' .AND. LEVEL_Q1 == 0 .AND. LEVEL_Q2 == 0) THEN
                LEVEL_P = LEVEL_P - 1
                IF (LEVEL_P == 0) THEN
                    JB = K - 1
                ENDIF
            ENDIF
            IF (TEMP_STATEMENT(K:K) == '"' .AND. LEVEL_Q1 == 0) LEVEL_Q2 = 1 - LEVEL_Q2
            IF (TEMP_STATEMENT(K:K) == "'" .AND. LEVEL_Q2 == 0) LEVEL_Q1 = 1 - LEVEL_Q1
            IF (TEMP_STATEMENT(K:K) == ',' .AND. LEVEL_Q1 == 0 .AND. LEVEL_Q2 == 0) THEN
                IF (LEVEL_P == 1) THEN
                    JB = K - 1
                ENDIF
            ENDIF
            IF (JB >= JA) THEN
                IF (R_ARGS(R_ARG_START(J)+NARG-1) <= 0) THEN
                    DO L = JA, JB-9
                       IF (TEMP_STATEMENT(L:L+6) == "TO_FM('") THEN
                           TEMP_STATEMENT(L:L+6) =  &
                           CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
                           DO L2 = L, JB
                              IF (TEMP_STATEMENT(L2:L2+1) == "')") THEN
                                  TEMP_STATEMENT(L2:L2+1) = CHAR(0)//CHAR(0)
                                  EXIT
                              ENDIF
                           ENDDO
                       ENDIF
                       IF (TEMP_STATEMENT(L:L+6) == "TO_IM('") THEN
                           TEMP_STATEMENT(L:L+6) =  &
                           CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
                           DO L2 = L, JB
                              IF (TEMP_STATEMENT(L2:L2+1) == "')") THEN
                                  TEMP_STATEMENT(L2:L2+1) = CHAR(0)//CHAR(0)
                                  EXIT
                              ENDIF
                           ENDDO
                       ENDIF
                    ENDDO
                ENDIF
                JA = K + 1
                NARG = NARG + 1
            ENDIF
            IF (LEVEL_P == 0) EXIT
         ENDDO

!             See if there are more references to this function in this statement.

         IF (JB > 1 .AND. JB < END_OF_STATEMENT) THEN
             KA = INDEX_ULC(TEMP_STATEMENT(JB:END_OF_STATEMENT),TRIM(R_NAME(J)))
             IF (KA > 0) THEN
                 KA = KA + JB - 1
                 GO TO 180
             ENDIF
         ENDIF
      ENDDO

!             Re-pack the statement if any TO_FM or TO_IM calls were removed.

      JA = 0
      DO K = 1, LAST_NONBLANK
         IF (TEMP_STATEMENT(K:K) /= CHAR(0)) THEN
             JA = JA + 1
             IF (K /= JA) TEMP_STATEMENT(JA:JA) = TEMP_STATEMENT(K:K)
         ENDIF
      ENDDO

      END_OF_STATEMENT = END_OF_STATEMENT - (LAST_NONBLANK - JA)
      LAST_NONBLANK = JA

      END SUBROUTINE CONVERT_CONSTANTS

      SUBROUTINE CONVERT_READS(LAST_NONBLANK,END_OF_STATEMENT)

!  TEMP_STATEMENT contains the statement.
!  Look for an executable statement where a read statement reads multiple precision variables,
!  and convert to the new form.
!
!  LAST_NONBLANK is the location of the last nonblank character in the statement,
!                including any trailing comments after !.
!  END_OF_STATEMENT is the location of the last nonblank character in the non-comment
!                   part of the statement.
!
!  The new read statement uses the same format as the original, and reads machine precision values.
!  Then new code converts to multiple precision.
!      READ (KR,*) (X(L),L=1,N)
!  becomes
!      ALLOCATE( C2FM_TEMP(1:N) )
!      READ (KR,*) (C2FM_TEMP(L),L=1,N)
!      DO L = 1, N
!         X(L) = C2FM_TEMP(L)
!      ENDDO
!      DEALLOCATE(C2FM_TEMP)
!  where X is type (fm) in the new version.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(63), PARAMETER :: NAME_CHARS =  &
        'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789'
      CHARACTER(32),  SAVE :: VAR_NAME
      CHARACTER(64),  SAVE :: READ_VAR_FM(100),READ_DP_VAR(100),READ_VAR
      CHARACTER(132), SAVE :: READ_SUB_FM1,READ_SUB_FM2
      INTEGER :: END_OF_STATEMENT,LAST_NONBLANK,ID1(6),ID2(6),J,J1,JA,JB,JC,JD,K,KA,KB,  &
                 K_VARIABLE,K_SECTION1,K_SECTION2,L,LAST_NB,P_LEVEL,Q1_LEVEL,Q2_LEVEL,INDEX_ULC,  &
                 N_VARS_FOUND,N_C2FM,BEGIN_READ_VAR_FM(100),END_READ_VAR_FM(100)
      LOGICAL :: IN_NAME,FM_VAR_FOUND

      N_AFTER_READ = 0
      TEMP_ST = TEMP_STATEMENT(1:END_OF_STATEMENT)
      KA = INDEX_ULC(TEMP_ST(1:END_OF_STATEMENT),'READ')
      IF (KA <= 0) RETURN

!             Scan backward from the keyword.  If there is a non-blank character before it,
!             make sure it is a digit (statement label) or right parentheses (end of logical if).

      DO J = KA-1, 1, -1
         IF (TEMP_ST(J:J) == ' ') CYCLE
         IF (INDEX('0123456789',TEMP_ST(J:J)) > 0) EXIT
         IF (TEMP_ST(J:J) == ')') EXIT
         RETURN
      ENDDO

!             Scan forward from the keyword.  The next field should be a format string or statement
!             label, followed by a comma ( '(e15.3,i5)', or 120, ), or a list with unit number,
!             format, ..., enclosed in parentheses ( (KW,120,end=555), or (6,"(e15.3,i5)") ).

      DO J = KA+4, END_OF_STATEMENT
         IF (TEMP_ST(J:J) == ' ') CYCLE
         IF (INDEX('0123456789',TEMP_ST(J:J)) > 0) THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == ' ') CYCLE
                IF (INDEX('0123456789',TEMP_ST(K:K)) > 0) CYCLE
                IF (TEMP_ST(K:K) == ',') THEN
                    KB = K + 1
                    GO TO 110
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == '"') THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == '"') THEN
                    DO L = K+1, END_OF_STATEMENT
                       IF (TEMP_ST(L:L) == ' ') CYCLE
                       IF (TEMP_ST(L:L) == ',') THEN
                           KB = L + 1
                           GO TO 110
                       ENDIF
                    ENDDO
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == "'") THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == "'") THEN
                    DO L = K+1, END_OF_STATEMENT
                       IF (TEMP_ST(L:L) == ' ') CYCLE
                       IF (TEMP_ST(L:L) == ',') THEN
                           KB = L + 1
                           GO TO 110
                       ENDIF
                    ENDDO
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == "(") THEN
             Q1_LEVEL = 0
             Q2_LEVEL = 0
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                IF (TEMP_ST(K:K) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                IF (TEMP_ST(K:K) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    KB = K + 1
                    GO TO 110
                ENDIF
             ENDDO
             RETURN
         ENDIF
      ENDDO
      RETURN

!             Now scan the i/o list starting at position KB and look for names of
!             multiple precision variables.

  110 JA = 0
      JB = 0
      N_AFTER_READ = 0
      N_VARS_FOUND = 0
      N_C2FM = 0
      READ_SUB_FM1 = ' '
      READ_SUB_FM2 = ' '
      LAST_NB = LAST_NONBLANK
      IN_NAME = .FALSE.
      FM_VAR_FOUND = .FALSE.
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      DO J = KB, END_OF_STATEMENT
         IF (TEMP_ST(J:J) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (TEMP_ST(J:J) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (Q1_LEVEL == 1 .OR. Q2_LEVEL == 1) CYCLE
         IF (.NOT. IN_NAME) THEN
             J1 = INDEX(NAME_CHARS,TEMP_ST(J:J))
             IF (J1 > 0) THEN
                 IN_NAME = .TRUE.
                 JA = J
                 IF (J == END_OF_STATEMENT) THEN
                     JB = J
                     VAR_NAME = ' '
                     CALL ALL_CAPS(TEMP_ST(JA:JB),JB-JA+1,VAR_NAME)
                     DO K = 1, N_VARS
                        IF (VAR_NAME == ROUTINE_VARIABLE(K)) THEN
                            FM_VAR_FOUND = .TRUE.
                            K_VARIABLE = K
                            EXIT
                        ENDIF
                     ENDDO
                 ENDIF
             ENDIF
         ELSE
             J1 = INDEX(NAME_CHARS,TEMP_ST(J:J))
             IF (J1 == 0 .OR. J == END_OF_STATEMENT) THEN
                 IN_NAME = .FALSE.
                 IF (J1 == 0) THEN
                     JB = J - 1
                 ELSE
                     JB = J
                 ENDIF
                 VAR_NAME = ' '
                 CALL ALL_CAPS(TEMP_ST(JA:JB),JB-JA+1,VAR_NAME)
                 DO K = 1, N_VARS
                    IF (VAR_NAME == ROUTINE_VARIABLE(K)) THEN
                        FM_VAR_FOUND = .TRUE.
                        K_VARIABLE = K
                        EXIT
                    ENDIF
                 ENDDO
             ENDIF
         ENDIF

         IF (FM_VAR_FOUND) THEN
             JC = JB

!             Check for a subscript after the name.

!             K_SECTION1 is positive if the "subscript" is really an array section.
!             Only 1 or 2 dimensional sections and only the simplest forms are handled.
!             For example, a( 2:n:2 ) or a( j , : ) are not recognized.

             P_LEVEL = 0
             K_SECTION1 = 0
             K_SECTION2 = 0
             ID1 = 0
             ID2 = 0
             DO K = JB+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == ' ') CYCLE
                IF (TEMP_ST(K:K) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                IF (TEMP_ST(K:K) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                IF (Q1_LEVEL == 1 .OR. Q2_LEVEL == 1) CYCLE
                IF (TEMP_ST(K:K) == ':' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
                    P_LEVEL == 1) THEN
                    IF (K_SECTION1 == 0) THEN
                        K_SECTION1 = K
                        ID1(4) = K - 1
                        ID1(5) = K + 1
                    ELSE IF (K_SECTION2 == 0) THEN
                        K_SECTION2 = K
                        ID1(4) = K - 1
                        ID1(5) = K + 1
                    ELSE
                        K_SECTION1 = 0
                    ENDIF
                    CYCLE
                ENDIF
                IF (TEMP_ST(K:K) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    P_LEVEL = P_LEVEL + 1
                    IF (P_LEVEL == 1) ID1(3) = K + 1
                    CYCLE
                ENDIF
                IF (P_LEVEL == 0) EXIT
                IF (TEMP_ST(K:K) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    P_LEVEL = P_LEVEL - 1
                    IF (P_LEVEL == 0) THEN
                        IF (K_SECTION2 > 0) THEN
                            ID2(6) = K - 1
                        ELSE IF (K_SECTION1 > 0) THEN
                            ID1(6) = K - 1
                        ENDIF
                        JC = K
                        EXIT
                    ENDIF
                ENDIF
             ENDDO

!             Look for an implied do after the subscript.
!             This will also handle two implied do loops after a 2-dimensional array.

             JD = JC
             IF (JC > JB .AND. K_SECTION1 == 0) THEN
                 DO K = JC+1, END_OF_STATEMENT
                    IF (TEMP_ST(K:K) == ' ') CYCLE
                    IF (TEMP_ST(K:K) == ',') THEN
                        CALL IMPLIED_DO(TEMP_ST,END_OF_STATEMENT,K,ID1,JD)
                        IF (JD > JC) THEN
                            DO L = JD+1, END_OF_STATEMENT
                               IF (TEMP_ST(L:L) == ' ') CYCLE
                               IF (TEMP_ST(L:L) == ',') THEN
                                   CALL IMPLIED_DO(TEMP_ST,END_OF_STATEMENT,L,ID2,JD)
                               ENDIF
                               EXIT
                            ENDDO
                        ENDIF
                    ENDIF
                    EXIT
                 ENDDO
             ENDIF

!             Replace the multiple precision variable name in the read list by a temporary double
!             precision variable, and generate new statements to convert the d.p. results to
!             multiple precision.

             IF (VARIABLE_TYPE(K_VARIABLE) == 1) THEN
                 N_VARS_FOUND = N_VARS_FOUND + 1
                 N_C2FM = N_C2FM + 1
                 IF (JC == JB) THEN
                     IF (READ_SUB_FM1(1:1) == ' ') THEN
                         READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                         BEGIN_READ_VAR_FM(N_C2FM) = JA
                         END_READ_VAR_FM(N_C2FM) = JC
                         READ_SUB_FM1(1:1) = '1'
                         READ_DP_VAR(N_C2FM) = 'C2FM_TEMP('//TRIM(READ_SUB_FM1)//')'
                     ELSE
                         READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                         BEGIN_READ_VAR_FM(N_C2FM) = JA
                         END_READ_VAR_FM(N_C2FM) = JC
                         READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+1'
                         READ_DP_VAR(N_C2FM) = 'C2FM_TEMP('//TRIM(READ_SUB_FM1)//')'
                     ENDIF
                 ELSE IF (K_SECTION1 <= 0) THEN
                     IF (ID1(1) <= 0) THEN
                         IF (READ_SUB_FM1(1:1) == ' ') THEN
                             READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                             BEGIN_READ_VAR_FM(N_C2FM) = JA
                             END_READ_VAR_FM(N_C2FM) = JC
                             READ_SUB_FM1(1:1) = '1'
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP('//TRIM(READ_SUB_FM1)//')'
                         ELSE
                             READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                             BEGIN_READ_VAR_FM(N_C2FM) = JA
                             END_READ_VAR_FM(N_C2FM) = JC
                             READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+1'
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP('//TRIM(READ_SUB_FM1)//')'
                         ENDIF
                     ELSE IF (ID2(1) <= 0) THEN
                         IF (READ_SUB_FM1(1:1) == ' ') THEN
                             READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                             BEGIN_READ_VAR_FM(N_C2FM) = JA
                             END_READ_VAR_FM(N_C2FM) = JC
                             IF (TEMP_ST(ID1(3):ID1(4)) == '1') THEN
                                 READ_SUB_FM1 = TEMP_ST(ID1(5):ID1(6))
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JC)
                             ELSE
                                 READ_SUB_FM1 = TEMP_ST(ID1(5):ID1(6))//'-'//  &
                                                     TEMP_ST(ID1(3):ID1(4))//'+1'
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JC)
                             ENDIF
                         ELSE
                             READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                             BEGIN_READ_VAR_FM(N_C2FM) = JA
                             END_READ_VAR_FM(N_C2FM) = JC
                             IF (TEMP_ST(ID1(3):ID1(4)) == '1') THEN
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JC-1)//'+'//  &
                                                       TRIM(READ_SUB_FM1)//')'
                                 READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+'//TEMP_ST(ID1(5):ID1(6))
                             ELSE
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JC-1)//'+'//  &
                                                       TRIM(READ_SUB_FM1)//')'
                                 READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+'//TEMP_ST(ID1(5):ID1(6))//  &
                                                '-'//TEMP_ST(ID1(3):ID1(4))//'+1'
                             ENDIF
                         ENDIF
                     ELSE
                         IF (READ_SUB_FM2(1:1) == ' ') THEN
                             READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                             BEGIN_READ_VAR_FM(N_C2FM) = JA
                             END_READ_VAR_FM(N_C2FM) = JC
                             IF (TEMP_ST(ID1(3):ID1(4)) == '1' .AND.  &
                                 TEMP_ST(ID2(3):ID2(4)) == '1') THEN
                                 READ_SUB_FM2 = TEMP_ST(ID1(5):ID1(6))//','//TEMP_ST(ID2(5):ID2(6))
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP2'//TEMP_ST(JB+1:JC)
                             ELSE
                                 READ_SUB_FM2 = TEMP_ST(ID1(5):ID1(6))//'-'//        &
                                                TEMP_ST(ID1(3):ID1(4))//'+1'//','//  &
                                                TEMP_ST(ID2(5):ID2(6))//'-'//        &
                                                TEMP_ST(ID2(3):ID2(4))//'+1'
                                 READ_DP_VAR(N_C2FM) = 'C2FM_TEMP2'//TEMP_ST(JB+1:JC)
                             ENDIF
                         ELSE
                             WRITE (*,*) ' '
                             WRITE (*,*) ' Convert2FM warning!    Only one 2-dimensional'//  &
                                         ' implied do can be handled in one read statement.'
                             WRITE (*,*) ' This statement seems to have more than one:'
                             WRITE (*,*) TRIM(TEMP_ST)
                             WRITE (*,*) ' '
                         ENDIF
                     ENDIF
                 ELSE IF (K_SECTION2 <= 0) THEN
                     IF (READ_SUB_FM1(1:1) == ' ') THEN
                         READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                         BEGIN_READ_VAR_FM(N_C2FM) = JA
                         END_READ_VAR_FM(N_C2FM) = JC
                         IF (TEMP_ST(ID1(3):ID1(4)) == '1') THEN
                             READ_SUB_FM1 = TEMP_ST(ID1(5):ID1(6))
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JC)
                         ELSE
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JB+1)//'1:'//  &
                                                   TEMP_ST(ID1(5):ID1(6))//'-'//  &
                                                   TEMP_ST(ID1(3):ID1(4))//'+1'//')'
                             READ_SUB_FM1 = TEMP_ST(ID1(5):ID1(6))//'-'//  &
                                            TEMP_ST(ID1(3):ID1(4))//'+1'
                         ENDIF
                     ELSE
                         READ_VAR_FM(N_C2FM) = TEMP_ST(JA:JC)
                         BEGIN_READ_VAR_FM(N_C2FM) = JA
                         END_READ_VAR_FM(N_C2FM) = JC
                         IF (TEMP_ST(ID1(3):ID1(4)) == '1') THEN
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JB+1)//  &
                                                   TRIM(READ_SUB_FM1)//'+'//TEMP_ST(JB+2:JC-1)//  &
                                                   '+'//TRIM(READ_SUB_FM1)//')'
                             READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+'//TEMP_ST(ID1(5):ID1(6))
                         ELSE
                             READ_DP_VAR(N_C2FM) = 'C2FM_TEMP'//TEMP_ST(JB+1:JB+1)//  &
                                                   TRIM(READ_SUB_FM1)//'+'//'1:'//  &
                                                   TEMP_ST(ID1(5):ID1(6))//'-'//  &
                                                   TEMP_ST(ID1(3):ID1(4))//'+1'//  &
                                                   '+'//TRIM(READ_SUB_FM1)//')'
                             READ_SUB_FM1 = TRIM(READ_SUB_FM1)//'+'//TEMP_ST(ID1(5):ID1(6))//  &
                                            '-'//TEMP_ST(ID1(3):ID1(4))//'+1'
                         ENDIF
                     ENDIF
                 ENDIF
                 CALL READ_VAR_SUBSCRIPT(READ_DP_VAR(N_C2FM))

!                Record the new conversion statements that will be written after this
!                read statement.

                 N_AFTER_READ = N_AFTER_READ + 1
                 IF (ID1(1) == 0) THEN
                     AFTER_READ(N_AFTER_READ) = '      '//TRIM(READ_VAR_FM(N_C2FM))//' = '//     &
                                                TRIM(READ_DP_VAR(N_C2FM))
                 ELSE IF (ID2(1) == 0) THEN
                     AFTER_READ(N_AFTER_READ) = '      '//'DO '//TEMP_ST(ID1(1):ID1(2))//        &
                                                ' = '//TEMP_ST(ID1(3):ID1(4))//', '//            &
                                                TEMP_ST(ID1(5):ID1(6))
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '         '//TRIM(READ_VAR_FM(N_C2FM))//' = '//  &
                                                TRIM(READ_DP_VAR(N_C2FM))
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '      ENDDO'
                 ELSE
                     AFTER_READ(N_AFTER_READ) = '      '//'DO '//TEMP_ST(ID1(1):ID1(2))//        &
                                                ' = '//TEMP_ST(ID1(3):ID1(4))//', '//            &
                                                TEMP_ST(ID1(5):ID1(6))
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '         '//'DO '//TEMP_ST(ID2(1):ID2(2))//     &
                                                ' = '//TEMP_ST(ID2(3):ID2(4))//', '//            &
                                                TEMP_ST(ID2(5):ID2(6))
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '            '//TRIM(READ_VAR_FM(N_C2FM))//      &
                                                ' = '//TRIM(READ_DP_VAR(N_C2FM))
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '         ENDDO'
                     N_AFTER_READ = N_AFTER_READ + 1
                     AFTER_READ(N_AFTER_READ) = '      ENDDO'
                 ENDIF
             ENDIF

             P_LEVEL = 0
             Q1_LEVEL = 0
             Q2_LEVEL = 0
             FM_VAR_FOUND = .FALSE.
         ENDIF
      ENDDO

!             If there are any new statements going after the read statement, put the corresponding
!             new statements allocating the temporary double precision arrays before the read.

      J1 = 0
      IF (N_AFTER_READ > 0) THEN
          IF (READ_SUB_FM1 /= ' ') THEN
              READ_VAR = 'C2FM_TEMP('//TRIM(READ_SUB_FM1)//')'
              CALL READ_VAR_SUBSCRIPT(READ_VAR)
              J1 = 1
              WRITE (23,"(A)") ' '
              WRITE (23,"(A)") '      ALLOCATE('//TRIM(READ_VAR)//')'
              N_AFTER_READ = N_AFTER_READ + 1
              AFTER_READ(N_AFTER_READ) = '      DEALLOCATE(C2FM_TEMP)'
          ENDIF
          IF (READ_SUB_FM2 /= ' ') THEN
              READ_VAR = 'C2FM_TEMP2('//TRIM(READ_SUB_FM2)//')'
              CALL READ_VAR_SUBSCRIPT(READ_VAR)
              IF (J1 == 0) WRITE (23,"(A)") ' '
              WRITE (23,"(A)") '      ALLOCATE('//TRIM(READ_VAR)//')'
              N_AFTER_READ = N_AFTER_READ + 1
              AFTER_READ(N_AFTER_READ) = '      DEALLOCATE(C2FM_TEMP2)'
          ENDIF
      ENDIF

!             Now that the entire statement has been scanned and the lists of old/new variables
!             built, change the variables in the read list.

      JA = 0
      KA = 1
      TEMP_STATEMENT = ' '
      J = 1
  120 IF (J > LAST_NONBLANK) GO TO 130
         JA = JA + 1
         IF (KA <= N_C2FM) THEN
             IF (J == BEGIN_READ_VAR_FM(KA)) THEN
                 JC = LEN_TRIM(READ_VAR_FM(KA))
                 JD = LEN_TRIM(READ_DP_VAR(KA))
                 TEMP_STATEMENT(JA:JA+JD-1) = TRIM(READ_DP_VAR(KA))
                 JA = JA + JD - 1
                 J = J + JC - 1
                 KA = KA + 1
                 LAST_NB = LAST_NB + JD - JC
             ELSE
                 TEMP_STATEMENT(JA:JA) = TEMP_ST(J:J)
             ENDIF
         ELSE
             TEMP_STATEMENT(JA:JA) = TEMP_ST(J:J)
         ENDIF
      J = J + 1
      IF (J <= LAST_NONBLANK) GO TO 120

  130 END_OF_STATEMENT = END_OF_STATEMENT + LAST_NB - LAST_NONBLANK
      LAST_NONBLANK = LAST_NB

      END SUBROUTINE CONVERT_READS

      SUBROUTINE CONVERT_TYPE_DECLARE(LAST_NONBLANK,END_OF_STATEMENT)

!  Scan a multiple precision type declaration statement, (in TEMP_STATEMENT), such as
!     TYPE (FM)        ,SAVE :: X(60),Y
!  and remove the SAVE attribute here.  Also remove extra blanks.  Later, a new SAVE statement
!  will be inserted, saving all the local multiple precision variables in this routine.
!     TYPE (FM) :: X(60),Y
!     ...
!     SAVE :: X,Y

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: LAST_NONBLANK,END_OF_STATEMENT
      INTEGER :: J,KA,KB,K1,INDEX_ULC

      KB = INDEX_ULC(TEMP_STATEMENT(1:END_OF_STATEMENT),'::')
      IF (KB < 1) RETURN
      KA = INDEX_ULC(TEMP_STATEMENT(1:KB),'SAVE')
      IF (KA > 0) THEN
          TEMP_STATEMENT(KA:KA+3) = '    '
          DO J = KA, 1, -1
             IF (TEMP_STATEMENT(J:J) == ',') THEN
                 TEMP_STATEMENT(J:J) = ' '
                 EXIT
             ENDIF
          ENDDO
      ENDIF
      KA = INDEX_ULC(TEMP_STATEMENT(1:KB),'TYPE (')
      IF (KA < 2) RETURN

      TEMP_ST(KA:END_OF_STATEMENT) = ' '
      TEMP_ST(1:KA) = TEMP_STATEMENT(1:KA)
      K1 = KA
      DO J = KA, KB
         IF (TEMP_STATEMENT(J-1:J) /= '  ') THEN
             TEMP_ST(K1:K1) = TEMP_STATEMENT(J:J)
             K1 = K1 + 1
         ENDIF
      ENDDO

      DO J = KB+1, LAST_NONBLANK
         TEMP_ST(K1:K1) = TEMP_STATEMENT(J:J)
         K1 = K1 + 1
      ENDDO
      K1 = K1 - 1

      TEMP_STATEMENT(1:LAST_NONBLANK) = ' '
      END_OF_STATEMENT = END_OF_STATEMENT - (LAST_NONBLANK - K1)
      LAST_NONBLANK = K1

      TEMP_STATEMENT(1:LAST_NONBLANK) = TEMP_ST(1:LAST_NONBLANK)

      END SUBROUTINE CONVERT_TYPE_DECLARE

      SUBROUTINE CONVERT_WRITES(LAST_NONBLANK,END_OF_STATEMENT)

!  TEMP_STATEMENT contains the statement.
!  Look for an executable statement where a write or print statement writes multiple precision
!  variables, and convert to the new form.
!
!  LAST_NONBLANK is the location of the last nonblank character in the statement,
!                including any trailing comments after !.
!  END_OF_STATEMENT is the location of the last nonblank character in the non-comment
!                   part of the statement.
!
!  The new write statement uses the same format as the original, and inserts a conversion
!  to machine precision (TO_DP in this example).
!      WRITE (KW,120) N,H,TOL
!  becomes
!      WRITE (KW,120) N,TO_DP(H),TO_DP(TOL)
!  where H and TOL are type (fm) in the new version.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(63), PARAMETER :: NAME_CHARS =  &
        'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789'
      CHARACTER(32), SAVE :: VAR_NAME
      INTEGER :: END_OF_STATEMENT,LAST_NONBLANK,J,J1,JA,JB,JC,K,KA,KB,K_VAR,L,  &
                 LAST_NB,P_LEVEL,Q1_LEVEL,Q2_LEVEL,INDEX_ULC,N_FCALLS,    &
                 FCALL_START(100),FCALL_END(100)
      LOGICAL :: IN_NAME,FM_VAR_FOUND

      TEMP_ST = TEMP_STATEMENT(1:END_OF_STATEMENT)
      KA = INDEX_ULC(TEMP_ST(1:END_OF_STATEMENT),'WRITE')
      KB = INDEX_ULC(TEMP_ST(1:END_OF_STATEMENT),'PRINT')
      IF (KA <= 0 .AND. KB <= 0) RETURN
      IF (KA > 0 .AND. KB > 0) THEN
          KA = MIN(KA,KB)
      ELSE IF (KB > 0) THEN
          KA = KB
      ENDIF

!             Scan backward from the keyword WRITE or PRINT.  If there is a non-blank character
!             before it, make sure it is a digit (statement label) or right parentheses (end of
!             logical if).

      DO J = KA-1, 1, -1
         IF (TEMP_ST(J:J) == ' ') CYCLE
         IF (INDEX('0123456789',TEMP_ST(J:J)) > 0) EXIT
         IF (TEMP_ST(J:J) == ')') EXIT
         RETURN
      ENDDO

!             Scan forward from the keyword.  The next field should be a format string or statement
!             label, followed by a comma ( '(e15.3,i5)', or 120, ), or a list with unit number,
!             format, ..., enclosed in parentheses ( (KW,120,end=555), or (6,"(e15.3,i5)") ).

      DO J = KA+5, END_OF_STATEMENT
         IF (TEMP_ST(J:J) == ' ') CYCLE
         IF (INDEX('0123456789',TEMP_ST(J:J)) > 0) THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == ' ') CYCLE
                IF (INDEX('0123456789',TEMP_ST(K:K)) > 0) CYCLE
                IF (TEMP_ST(K:K) == ',') THEN
                    KB = K + 1
                    GO TO 110
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == '"') THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == '"') THEN
                    DO L = K+1, END_OF_STATEMENT
                       IF (TEMP_ST(L:L) == ' ') CYCLE
                       IF (TEMP_ST(L:L) == ',') THEN
                           KB = L + 1
                           GO TO 110
                       ENDIF
                    ENDDO
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == "'") THEN
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == "'") THEN
                    DO L = K+1, END_OF_STATEMENT
                       IF (TEMP_ST(L:L) == ' ') CYCLE
                       IF (TEMP_ST(L:L) == ',') THEN
                           KB = L + 1
                           GO TO 110
                       ENDIF
                    ENDDO
                ENDIF
             ENDDO
             RETURN
         ENDIF
         IF (TEMP_ST(J:J) == "(") THEN
             Q1_LEVEL = 0
             Q2_LEVEL = 0
             DO K = J+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                IF (TEMP_ST(K:K) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                IF (TEMP_ST(K:K) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    KB = K + 1
                    GO TO 110
                ENDIF
             ENDDO
             RETURN
         ENDIF
      ENDDO
      RETURN

!             Look for multiprecision function calls in the write list.

  110 N_FCALLS = 0
      DO J = 1, N_ROUTINES
         JA = INDEX_ULC(TEMP_ST(1:END_OF_STATEMENT),TRIM(R_NAME(J)))
         IF (JA <= 0) CYCLE
         DO K = JA, END_OF_STATEMENT
            L = LEN_TRIM(R_NAME(J))
            IF (TEMP_ST(K:K+L-1) /= R_NAME(J)(1:L)) CYCLE
            IF (K > 1) THEN
                IF (INDEX(NAME_CHARS,TEMP_ST(K-1:K-1)) > 0) CYCLE
            ENDIF
            IF (K+L <= END_OF_STATEMENT) THEN
                IF (INDEX(NAME_CHARS,TEMP_ST(K+L:K+L)) > 0) CYCLE
            ENDIF
            N_FCALLS = N_FCALLS + 1
            FCALL_START(N_FCALLS) = K
            CALL SCAN_LIST(TEMP_ST,K,END_OF_STATEMENT,FCALL_END(N_FCALLS))
         ENDDO
      ENDDO

!             Now scan the i/o list starting at position KB and look for names of
!             multiple precision variables.
!             Scan backward so that the next match in TEMP_ST gives the same position in
!             TEMP_STATEMENT.  TEMP_STATEMENT is the one to be changed.

      JA = 0
      JB = 0
      LAST_NB = LAST_NONBLANK
      IN_NAME = .FALSE.
      FM_VAR_FOUND = .FALSE.
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      DO J = END_OF_STATEMENT, KB, -1
         IF (TEMP_ST(J:J) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (TEMP_ST(J:J) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (Q1_LEVEL == 1 .OR. Q2_LEVEL == 1) CYCLE
         IF (.NOT. IN_NAME) THEN
             J1 = INDEX(NAME_CHARS,TEMP_ST(J:J))
             IF (J1 > 0) THEN
                 IN_NAME = .TRUE.
                 JB = J
             ENDIF
         ELSE
             J1 = INDEX(NAME_CHARS,TEMP_ST(J:J))
             IF (J1 == 0) THEN
                 IN_NAME = .FALSE.
                 JA = J + 1
                 VAR_NAME = ' '
                 CALL ALL_CAPS(TEMP_ST(JA:JB),JB-JA+1,VAR_NAME)
                 DO K = 1, N_VARS
                    IF (VAR_NAME == ROUTINE_VARIABLE(K)) THEN
                        FM_VAR_FOUND = .TRUE.
                        K_VAR = K
                        EXIT
                    ENDIF
                 ENDDO
             ENDIF
         ENDIF

         IF (N_FCALLS > 0) THEN
             JC = 0
             DO K = 1, N_FCALLS
                IF (J >= FCALL_START(K) .AND. J <= FCALL_END(K)) JC = 1
             ENDDO
             IF (JC == 1) CYCLE
         ENDIF

         IF (FM_VAR_FOUND) THEN
             JC = JB

!             Check for a subscript or argument list after the name

             P_LEVEL = 0
             DO K = JB+1, END_OF_STATEMENT
                IF (TEMP_ST(K:K) == ' ') CYCLE
                IF (TEMP_ST(K:K) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                IF (TEMP_ST(K:K) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                IF (Q1_LEVEL == 1 .OR. Q2_LEVEL == 1) CYCLE
                IF (TEMP_ST(K:K) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    P_LEVEL = P_LEVEL + 1
                    CYCLE
                ENDIF
                IF (P_LEVEL == 0) EXIT
                IF (TEMP_ST(K:K) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                    P_LEVEL = P_LEVEL - 1
                    IF (P_LEVEL == 0) THEN
                        JC = K
                        EXIT
                    ENDIF
                ENDIF
             ENDDO

!             Enclose the item in TO_DP, TO_INT. or TO_DPZ.

             IF (VARIABLE_TYPE(K_VAR) == 1) THEN
                 LAST_NB = LAST_NB + 7
                 DO K = LAST_NB, JC+8, -1
                    TEMP_STATEMENT(K:K) = TEMP_STATEMENT(K-7:K-7)
                 ENDDO
                 TEMP_STATEMENT(JC+7:JC+7) = ')'
                 DO K = JC, JA, -1
                    TEMP_STATEMENT(K+6:K+6) = TEMP_STATEMENT(K:K)
                 ENDDO
                 TEMP_STATEMENT(JA:JA+5) = 'TO_DP('
             ELSE IF (VARIABLE_TYPE(K_VAR) == 2) THEN
                 LAST_NB = LAST_NB + 8
                 DO K = LAST_NB, JC+9, -1
                    TEMP_STATEMENT(K:K) = TEMP_STATEMENT(K-8:K-8)
                 ENDDO
                 TEMP_STATEMENT(JC+8:JC+8) = ')'
                 DO K = JC, JA, -1
                    TEMP_STATEMENT(K+7:K+7) = TEMP_STATEMENT(K:K)
                 ENDDO
                 TEMP_STATEMENT(JA:JA+6) = 'TO_DPZ('
             ELSE IF (VARIABLE_TYPE(K_VAR) == 3) THEN
                 LAST_NB = LAST_NB + 8
                 DO K = LAST_NB, JC+9, -1
                    TEMP_STATEMENT(K:K) = TEMP_STATEMENT(K-8:K-8)
                 ENDDO
                 TEMP_STATEMENT(JC+8:JC+8) = ')'
                 DO K = JC, JA, -1
                    TEMP_STATEMENT(K+7:K+7) = TEMP_STATEMENT(K:K)
                 ENDDO
                 TEMP_STATEMENT(JA:JA+6) = 'TO_INT('
             ENDIF
             P_LEVEL = 0
             Q1_LEVEL = 0
             Q2_LEVEL = 0
             FM_VAR_FOUND = .FALSE.
         ENDIF
      ENDDO

      END_OF_STATEMENT = END_OF_STATEMENT + LAST_NB - LAST_NONBLANK
      LAST_NONBLANK = LAST_NB

!             Look for constants to change inside multiprecision function calls in the write list.
!             The statement may have changed, so get FCALL_START and FCALL_END again.

      N_FCALLS = 0
      DO J = 1, N_ROUTINES
         JA = INDEX_ULC(TEMP_STATEMENT(1:END_OF_STATEMENT),TRIM(R_NAME(J)))
         IF (JA <= 0) CYCLE
         DO K = JA, END_OF_STATEMENT
            L = LEN_TRIM(R_NAME(J))
            IF (TEMP_STATEMENT(K:K+L-1) /= R_NAME(J)(1:L)) CYCLE
            IF (K > 1) THEN
                IF (INDEX(NAME_CHARS,TEMP_STATEMENT(K-1:K-1)) > 0) CYCLE
            ENDIF
            IF (K+L <= END_OF_STATEMENT) THEN
                IF (INDEX(NAME_CHARS,TEMP_STATEMENT(K+L:K+L)) > 0) CYCLE
            ENDIF
            N_FCALLS = N_FCALLS + 1
            FCALL_START(N_FCALLS) = K
            CALL SCAN_LIST(TEMP_STATEMENT,K,END_OF_STATEMENT,FCALL_END(N_FCALLS))
         ENDDO
      ENDDO
      DO J = N_FCALLS, 1, -1
         CALL CONVERT_CONSTANTS(1,LAST_NONBLANK,END_OF_STATEMENT,FCALL_START(J),FCALL_END(J))
      ENDDO

      END SUBROUTINE CONVERT_WRITES

      SUBROUTINE EXPAND_USE_LIST(J,K)

!  Routine J uses module K.  Check the list of used modules for module K
!  to see if any modules used there need to be added to the R_ALL_USED list.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: J,K,L,M

      IF (R_N_ALL_USED(K) <= 0) RETURN
      DO L = 1, R_N_ALL_USED(K)
         DO M = 1, R_N_ALL_USED(J)
            IF (R_ALL_USED(J,M) == R_ALL_USED(K,L)) EXIT
            IF (M == R_N_ALL_USED(J)) THEN
                R_N_ALL_USED(J) = R_N_ALL_USED(J) + 1
                R_ALL_USED(J,R_N_ALL_USED(J)) = R_ALL_USED(K,L)
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE EXPAND_USE_LIST

      SUBROUTINE FIND_NAME(NAME,LIST,N_ROUTINES,J2)

!  Look up NAME on LIST( 1 : N_ROUTINES ).
!  J2 is returned positive if LIST(J2) is NAME.

      IMPLICIT NONE
      INTEGER :: J,J2,N_ROUTINES
      CHARACTER(32) :: NAME,LIST(N_ROUTINES)

      J2 = 0
      DO J = 1, N_ROUTINES
         IF (LIST(J) == NAME) THEN
             J2 = J
             RETURN
         ENDIF
      ENDDO

      END SUBROUTINE FIND_NAME

      SUBROUTINE FIRST_EXECUTABLE(LAST,J1,J2)

!  The current routine is in LINE(1), ..., LINE(LAST).
!  J1 is returned as the line after the end of a PROGRAM, SUBROUTINE, FUNCTION, or MODULE statement.
!  J2 is returned as the first line of the first executable statement in the current routine.
!     Some new statements may need to be inserted there.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132), SAVE :: TEMP,PREVIOUS_LINE
      INTEGER :: FIRST,FIRST_NONBLANK,J,J1,J2,K,K1,KA,KL(4),LAST
      INTEGER, EXTERNAL :: INDEX_ULC

      J1 = 1
      J2 = LAST + 1
      DO J = 1, LAST
         TEMP = LINE(J)
         K1 = FIRST_NONBLANK(TEMP)
         IF (K1 < 1) CYCLE

!             Check for the first statement of the routine.

         KL(1) = INDEX_ULC(TEMP,'PROGRAM')
         KL(2) = INDEX_ULC(TEMP,'MODULE')
         KL(3) = INDEX_ULC(TEMP,'SUBROUTINE')
         KL(4) = INDEX_ULC(TEMP,'FUNCTION')
         K = 0
         DO K1 = 1, 4
            IF (K == 0 .AND. KL(K1) > 0) K = KL(K1)
            IF (K > 0  .AND. KL(K1) > 0) K = MIN(K,KL(K1))
         ENDDO
         IF (K > 0) CALL FIRST_STATEMENT(TEMP,FIRST,1,K)
         IF (FIRST > 0) THEN
             J1 = J + 1

!             Check for continuation lines.

             DO K1 = J, LAST
                TEMP = LINE(K1)
                KA = 0
                DO K = 1, 132
                   IF (TEMP(K:K) /= ' ') THEN
                       IF (TEMP(K:K) /= '!') KA = 1
                       EXIT
                   ENDIF
                ENDDO
                IF (KA == 0) CYCLE
                KA = 1
                DO K = 1, 132
                   IF (TEMP(K:K) == '!') EXIT
                   IF (TEMP(K:K) /= ' ') KA = K
                ENDDO
                IF (TEMP(KA:KA) /= '&') THEN
                    J1 = K1 + 1
                    GO TO 110
                ENDIF
             ENDDO
         ENDIF
      ENDDO

!             Scan past any declaration statements to find the first executable statement.

  110 TEMP = ' '
      DO J = J1, LAST
         K1 = FIRST_NONBLANK(LINE(J))
         IF (K1 < 1) CYCLE
         IF (LINE(J)(K1:K1) == '!') CYCLE
         PREVIOUS_LINE = TEMP
         TEMP = LINE(J)
         IF (INDEX_ULC(TEMP,'ALLOCATABLE') == K1) CYCLE
         IF (INDEX_ULC(TEMP,'CHARACTER'  ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'COMPLEX'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'DATA'       ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'DIMENSION'  ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'DOUBLE'     ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'EXTERNAL'   ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'FUNCTION'   ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'IMPLICIT'   ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'INTEGER'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'INTENT'     ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'LOGICAL'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'MODULE'     ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'NAMELIST'   ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'OPTIONAL'   ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'PARAMETER'  ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'POINTER'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'PROGRAM'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'PRIVATE'    ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'PUBLIC'     ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'REAL'       ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'SAVE'       ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'SUBROUTINE' ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'TYPE'       ) == K1) CYCLE
         IF (INDEX_ULC(TEMP,'USE'        ) == K1) CYCLE
         IF (J > J1) THEN
             IF (INDEX(PREVIOUS_LINE,'&') > 0) THEN
                 K1 = 0
                 KA = INDEX(PREVIOUS_LINE,'!')
                 IF (KA == 0) KA = 133
                 DO K = KA-1, 1, -1
                    IF (PREVIOUS_LINE(K:K) == ' ') CYCLE
                    IF (PREVIOUS_LINE(K:K) == '&') THEN
                        K1 = 1
                        EXIT
                    ENDIF
                 ENDDO
                 IF (K1 == 1) CYCLE
             ENDIF
         ENDIF
         J2 = J
         RETURN
      ENDDO

      END SUBROUTINE FIRST_EXECUTABLE

      FUNCTION FIRST_NONBLANK(LINE)

!  Find the first nonblank character on the line.
!  Digits (statement labels) are skipped, along with blanks.

      IMPLICIT NONE
      CHARACTER(132) :: LINE
      INTEGER :: FIRST_NONBLANK,J

      FIRST_NONBLANK = 0
      DO J = 1, 132
         IF (INDEX(' 0123456789',LINE(J:J)) <= 0) THEN
             FIRST_NONBLANK = J
             RETURN
         ENDIF
      ENDDO

      END FUNCTION FIRST_NONBLANK

      SUBROUTINE FIRST_STATEMENT(LINE1,FIRST,J,K)

!  Check to see if LINE1 is the start of the PROGRAM, SUBROUTINE, FUNCTION, or MODULE subprogram.
!  One of those four keywords starts in position K of LINE1.
!  FIRST is returned zero if LINE1 is not one of those statements, J if it is.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132) :: LINE1
      INTEGER :: FIRST,J,J1,J2,J3,K,L,JL1,INDEX_ULC

      FIRST = 0
      IF (INDEX(LINE1(1:K),'!') > 0) RETURN

!                    Check for a MODULE PROCEDURE statement.

      J1 = INDEX_ULC(LINE1(1:132),'MODULE')
      J2 = INDEX_ULC(LINE1(1:132),'PROCEDURE')
      J3 = INDEX(LINE1(1:132),'!')
      IF (J3 == 0) J3 = 132
      IF (J1 > 0 .AND. J2 > 0 .AND. J1 < J2 .AND. J1 < J3 .AND. J2 < J3) RETURN

      IF (K >= 1) THEN
          DO L = 1, K
             IF (LINE1(L:L) /= ' ') THEN
                 IF (INDEX_ULC(LINE1(L:L+6),'PROGRAM') > 0 .AND. LINE1(L+7:L+7) == ' ') THEN
                     FIRST = J
                     RETURN
                 ENDIF
                 IF (INDEX_ULC(LINE1(L:L+5),'MODULE') > 0 .AND. LINE1(L+6:L+6) == ' ') THEN
                     FIRST = J
                     RETURN
                 ENDIF
                 IF (INDEX_ULC(LINE1(L:L+9),'SUBROUTINE') > 0 .AND. LINE1(L+10:L+10) == ' ') THEN
                     FIRST = J
                     RETURN
                 ENDIF
                 IF (INDEX_ULC(LINE1(L:L+7),'FUNCTION') > 0 .AND. LINE1(L+8:L+8) == ' ') THEN
                     FIRST = J
                     RETURN
                 ENDIF
                 J1 = INDEX_ULC(LINE1(1:132),'RECURSIVE')
                 IF (J1 == K) THEN
                     J2 = INDEX_ULC(LINE1(1:132),'SUBROUTINE')
                     IF (J2 > K .AND. J2 < J3 .AND. LINE1(J1+9:J1+9) == ' ' .AND.  &
                         LINE1(J2-1:J2-1) == ' ' .AND. LINE1(J2+10:J2+10) == ' ') THEN
                         FIRST = J
                         RETURN
                     ENDIF
                     J2 = INDEX_ULC(LINE1(1:132),'FUNCTION')
                     IF (J2 > K .AND. J2 < J3 .AND. LINE1(J1+9:J1+9) == ' ' .AND.  &
                         LINE1(J2-1:J2-1) == ' ' .AND. LINE1(J2+8:J2+8) == ' ') THEN
                         FIRST = J
                         RETURN
                     ENDIF
                 ENDIF

                 JL1 = 0
                 J1 = INDEX_ULC(LINE1(1:132),'INTEGER')
                 IF (J1 == K) JL1 = 7
                 J1 = INDEX_ULC(LINE1(1:132),'REAL')
                 IF (J1 == K) JL1 = 4
                 J1 = INDEX_ULC(LINE1(1:132),'COMPLEX')
                 IF (J1 == K) JL1 = 7
                 J1 = INDEX_ULC(LINE1(1:132),'CHARACTER')
                 IF (J1 == K) JL1 = 9
                 J1 = INDEX_ULC(LINE1(1:132),'DOUBLE')
                 IF (J1 == K) JL1 = 6
                 J1 = INDEX_ULC(LINE1(1:132),'LOGICAL')
                 IF (J1 == K) JL1 = 7

                 IF (JL1 > 0) THEN
                     J1 = K
                     J2 = INDEX_ULC(LINE1(1:132),'SUBROUTINE')
                     IF (J2 > K .AND. J2 < J3 .AND. LINE1(J1+JL1:J1+JL1) == ' ' .AND.  &
                         LINE1(J2-1:J2-1) == ' ' .AND. LINE1(J2+10:J2+10) == ' ') THEN
                         FIRST = J
                         RETURN
                     ENDIF
                     J2 = INDEX_ULC(LINE1(1:132),'FUNCTION')
                     IF (J2 > K .AND. J2 < J3 .AND. LINE1(J1+JL1:J1+JL1) == ' ' .AND.  &
                         LINE1(J2-1:J2-1) == ' ' .AND. LINE1(J2+8:J2+8) == ' ') THEN
                         FIRST = J
                         RETURN
                     ENDIF
                 ENDIF
                 RETURN
             ENDIF
          ENDDO
      ENDIF

      END SUBROUTINE FIRST_STATEMENT

      SUBROUTINE GET_ARGUMENTS

!  Make a first pass over the program to find out which arguments to each subprogram
!  will be multiple precision.

!  ROUTINE_NAME         is the name of the current routine.
!  ROUTINE_ARGUMENT(K)  is the variable name of the Kth argument in the current routine.
!  ARGUMENT_TYPE(K) = 0 if that variable will not become multiple precision
!                   = 1 for TYPE (FM)
!                   = 2 for TYPE (ZM)
!                   = 3 for TYPE (IM)
!  N_ROUTINES           is the number of routines in the program.
!  R_NAME(J)            is the routine name of the Jth subprogram.
!  R_TYPE(J)            is the routine type of the Jth subprogram.
!  R_N_ARGS(J)          is the number of arguments for the Jth routine.
!  R_ARG_START(J)       gives the starting location of the list of argument types for routine J.
!  R_ARGS( R_ARG_START(J) : R_ARG_START(J) + R_N_ARGS(J) - 1 )
!                       gives the argument type of each input argument to the Jth routine,
!                       using the same codes as ARGUMENT_TYPE.
!  R_N_USED(J)          is the number of USE statements in the Jth routine.
!  R_USED_START(J)      gives the starting location of the list of used modules for routine J.
!  R_USED( R_USED_START(J) : R_USED_START(J) + R_N_USED(J) - 1 )
!                       gives the module names used in the Jth routine.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(32)   :: TEMP
      CHARACTER(132)  :: T_LINE
      INTEGER :: FIRST,FIRST_NONBLANK,J,K,K1,KA,KB,KL(4),L,L1,L2,LAST,TOTAL_VARS,  &
                 Q1_LEVEL,Q2_LEVEL,P_LEVEL,INDEX_ULC
      LOGICAL :: FILE_ENDED

      DO J = 1, 1000

!             Read the next routine.

         CALL READ_ROUTINE(22,FIRST,LAST,FILE_ENDED)

         IF (FILE_ENDED .AND. (FIRST <= 0 .OR. LAST <= 0)) RETURN
         N_ROUTINES = J
         R_N_USED(J) = 0
         IF (J <= 1) THEN
             R_USED_START(J) = 1
         ELSE
             R_USED_START(J) = R_USED_START(J-1) + R_N_USED(J-1)
         ENDIF

!             Scan the arguments on the first line.

!             Q1_LEVEL or Q2_LEVEL is positive when the scan in within a character string.
!             P_LEVEL is positive when the scan in within a pair of parentheses.

         Q1_LEVEL = 0
         Q2_LEVEL = 0
         P_LEVEL = 0
         N_ARGS = 0
         ROUTINE_ARGUMENT(1) = ' '
         ARGUMENT_TYPE(1:100) = -1

         KL(1) = INDEX_ULC(LINE(FIRST),'PROGRAM')
         KL(2) = INDEX_ULC(LINE(FIRST),'MODULE')
         KL(3) = INDEX_ULC(LINE(FIRST),'SUBROUTINE')
         KL(4) = INDEX_ULC(LINE(FIRST),'FUNCTION')
         K = 0
         DO K1 = 1, 4
            IF (K == 0 .AND. KL(K1) > 0) K = KL(K1)
            IF (K > 0  .AND. KL(K1) > 0) K = MIN(K,KL(K1))
         ENDDO
         KA = 0
         KB = K
         ROUTINE_NAME = ' '
         DO L = K, 132
            IF (LINE(FIRST)(L:L) == ' ' .AND. KA == 0) KA = 1
            IF (LINE(FIRST)(L:L) /= ' ' .AND. KA == 1) KA = L
            IF ((LINE(FIRST)(L:L) == ' ' .OR. LINE(FIRST)(L:L) == '(') .AND. KA > 1) THEN
                KB = L - 1
                CALL ALL_CAPS(LINE(FIRST)(KA:KB),KB-KA+1,ROUTINE_NAME)
                EXIT
             ENDIF
         ENDDO

         IF (INDEX_ULC(LINE(FIRST),'PROGRAM') > 0) GO TO 120
         IF (INDEX_ULC(LINE(FIRST),'MODULE') > 0) GO TO 120
         IF (K <= 0) GO TO 130

  110    DO L = KB, 132
            IF (LINE(FIRST)(L:L) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                P_LEVEL = P_LEVEL + 1
                IF (P_LEVEL == 1) KA = L + 1
            ENDIF
            IF (LINE(FIRST)(L:L) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                P_LEVEL = P_LEVEL - 1
                IF (P_LEVEL == 0) THEN
                    KB = L - 1
                    N_ARGS = N_ARGS + 1
                    ROUTINE_ARGUMENT(N_ARGS) = LINE(FIRST)(KA:KB)
                    EXIT
                ENDIF
            ENDIF
            IF (LINE(FIRST)(L:L) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
            IF (LINE(FIRST)(L:L) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
            IF (LINE(FIRST)(L:L) == ',' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
                P_LEVEL == 1) THEN
                KB = L - 1
                N_ARGS = N_ARGS + 1
                ROUTINE_ARGUMENT(N_ARGS) = LINE(FIRST)(KA:KB)
                KA = L + 1
            ENDIF
            IF (LINE(FIRST)(L:L) == '&' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
                P_LEVEL == 1) THEN
                FIRST = FIRST + 1
                KA = 1
                KB = 1
                GO TO 110
            ENDIF
         ENDDO

!             Remove leading blanks from the routine arguments.

         DO L = 1, N_ARGS
            DO KA = 1, 32
               IF (ROUTINE_ARGUMENT(L)(KA:KA) /= ' ') THEN
                   IF (KA > 1) THEN
                       TEMP = ROUTINE_ARGUMENT(L)
                       ROUTINE_ARGUMENT(L) = ' '
                       ROUTINE_ARGUMENT(L) = TEMP(KA:32)
                       EXIT
                   ELSE
                       EXIT
                   ENDIF
               ENDIF
            ENDDO
         ENDDO

!             Write the argument information to file C2FM.VARS for ordinary routines, or
!             to file C2FM.MODS for modules.

  120 IF (ROUTINE_TYPE == 4) THEN
          WRITE (32,"(//A,I4,A,A/)") ' Routine ',J,':   Module     ',TRIM(ROUTINE_NAME)
          WRITE (33,"(//A,I4,A,A/)") ' Routine ',J,':   Module     ',TRIM(ROUTINE_NAME)
      ELSE IF (ROUTINE_TYPE == 3) THEN
          WRITE (32,"(//A,I4,A,A/)") ' Routine ',J,':   Function   ',TRIM(ROUTINE_NAME)
      ELSE IF (ROUTINE_TYPE == 2) THEN
          WRITE (32,"(//A,I4,A,A/)") ' Routine ',J,':   Subroutine ',TRIM(ROUTINE_NAME)
      ELSE IF (ROUTINE_TYPE == 1) THEN
          WRITE (32,"(//A,I4,A,A/)") ' Routine ',J,':   Program    ',TRIM(ROUTINE_NAME)
      ENDIF

!             The SET_TYPE routine writes the information about all multiple precision variables
!             used in each routine.

!             Check to see which of these arguments are multiple-precision, and set ARGUMENT_TYPE.

         DO L = 1, N_ARGS
            ARGUMENT_TYPE(L) = 0
         ENDDO

         TOTAL_VARS = 0
         DO L = FIRST, LAST
            IF (INDEX_ULC(LINE(L),'REAL (KIND(3.1D1))') > 0) THEN
                CALL SET_TYPE(1,L,TOTAL_VARS)
                K = INDEX_ULC(LINE(L),'REAL (KIND(3.1D1))')
                LINE(L)(K:K+17) = 'TYPE (FM)         '
                CYCLE
            ENDIF
            IF (INDEX_ULC(LINE(L),'COMPLEX (KIND(3.1D1))') > 0) THEN
                CALL SET_TYPE(2,L,TOTAL_VARS)
                K = INDEX_ULC(LINE(L),'COMPLEX (KIND(3.1D1))')
                LINE(L)(K:K+20) = 'TYPE (ZM)            '
                CYCLE
            ENDIF
            IF (INDEX_ULC(LINE(L),'INTEGER (KIND(31))') > 0) THEN
                CALL SET_TYPE(3,L,TOTAL_VARS)
                K = INDEX_ULC(LINE(L),'INTEGER (KIND(31))')
                LINE(L)(K:K+17) = 'TYPE (IM)         '
                CYCLE
            ENDIF
         ENDDO

!             Write the argument information to file C2FM.ARGS, and build the database for
!             routine argument information.

  130    IF (ROUTINE_TYPE == 4) THEN
             WRITE (31,"(//A,I4,A,A/)") ' Routine ',J,':   Module     ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 3) THEN
             WRITE (31,"(//A,I4,A,A/)") ' Routine ',J,':   Function   ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 2) THEN
             WRITE (31,"(//A,I4,A,A/)") ' Routine ',J,':   Subroutine ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 1) THEN
             WRITE (31,"(//A,I4,A,A/)") ' Routine ',J,':   Program    ',TRIM(ROUTINE_NAME)
         ENDIF
         R_NAME(N_ROUTINES) = ROUTINE_NAME
         R_TYPE(N_ROUTINES) = ROUTINE_TYPE
         R_N_ARGS(N_ROUTINES) = N_ARGS
         IF (N_ROUTINES == 1) THEN
             R_ARG_START(N_ROUTINES) = 1
         ELSE
             R_ARG_START(N_ROUTINES) = R_ARG_START(N_ROUTINES-1) + R_N_ARGS(N_ROUTINES-1)
         ENDIF
         DO K = 1, N_ARGS
            IF (ARGUMENT_TYPE(K) == 0) THEN
                WRITE (31,"(A,I4,A,A/)") ' Argument ',K,'   Type: ----  ,  Name = ',  &
                                         TRIM(ROUTINE_ARGUMENT(K))
            ELSE IF (ARGUMENT_TYPE(K) == 1) THEN
                WRITE (31,"(A,I4,A,A/)") ' Argument ',K,'   Type: (FM)  ,  Name = ',  &
                                         TRIM(ROUTINE_ARGUMENT(K))
            ELSE IF (ARGUMENT_TYPE(K) == 2) THEN
                WRITE (31,"(A,I4,A,A/)") ' Argument ',K,'   Type: (ZM)  ,  Name = ',  &
                                         TRIM(ROUTINE_ARGUMENT(K))
            ELSE IF (ARGUMENT_TYPE(K) == 3) THEN
                WRITE (31,"(A,I4,A,A/)") ' Argument ',K,'   Type: (IM)  ,  Name = ',  &
                                         TRIM(ROUTINE_ARGUMENT(K))
            ENDIF
            R_ARGS(R_ARG_START(N_ROUTINES)+K-1) = ARGUMENT_TYPE(K)
         ENDDO

!             Check for USE statements in this routine.

         CALL FIRST_EXECUTABLE(LAST,KA,KB)
         DO K = KA, KB
            CALL ALL_CAPS(LINE(K),132,T_LINE)
            K1 = FIRST_NONBLANK(T_LINE)
            IF (K1 < 1) CYCLE
            IF (INDEX_ULC(T_LINE,'USE ') == K1) THEN
                L1 = 0
                L2 = 0
                DO L = K1+4, 132
                   IF (T_LINE(L:L) /= ' ') THEN
                       L1 = L
                       EXIT
                   ENDIF
                ENDDO
                IF (L1 > 0) THEN
                    DO L = L1, 132
                       IF (T_LINE(L:L) == ' ' .OR. T_LINE(L:L) == ',') THEN
                           L2 = L
                           EXIT
                       ENDIF
                    ENDDO
                ENDIF
                IF (L1 > 0 .AND. L2 >= L1) THEN
                    R_N_USED(J) = R_N_USED(J) + 1
                    R_USED(R_USED_START(J) + R_N_USED(J) - 1) = T_LINE(L1:L2)
                ENDIF
            ENDIF
         ENDDO

!             Write the preliminary modifications to file C2FM.PASS1.

         DO K = 1, LAST
            WRITE (34,"(A)") TRIM(LINE(K))
         ENDDO

         IF (FILE_ENDED) RETURN

      ENDDO

      END SUBROUTINE GET_ARGUMENTS

      FUNCTION GET_SIZE(STRING,KA,KB)

!  STRING(KA:KB) is a dimension declaration.  Find the total size of the array.

      IMPLICIT NONE
      CHARACTER(132) :: STRING
      INTEGER :: GET_SIZE,KA,KB,J,JSTART,JEND,KVAL

      GET_SIZE = 1
      JSTART = 1
      JEND = 0
      DO J = KA+1, KB
         IF (STRING(J:J) == ' ') CYCLE
         KVAL = INDEX('0123456789',STRING(J:J))
         IF (KVAL > 0) JEND = 10*JEND + (KVAL-1)
         IF (STRING(J:J) == ':') THEN
             JSTART = JEND
             JEND = 0
             CYCLE
         ENDIF
         IF (STRING(J:J) == ',') THEN
             GET_SIZE = GET_SIZE * (JEND-JSTART+1)
             JSTART = 1
             JEND = 0
             CYCLE
         ENDIF
         IF (STRING(J:J) == ')') THEN
             GET_SIZE = GET_SIZE * (JEND-JSTART+1)
             RETURN
         ENDIF
         IF (KVAL <= 0) THEN
             GET_SIZE = 0
             RETURN
         ENDIF
      ENDDO

      END FUNCTION GET_SIZE

      SUBROUTINE GET_STATEMENT(J,JB,START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT)

!  LINE(J) contains the start of the next statement.
!  Build the entire statement, including any continuation lines, in TEMP_STATEMENT.
!  JB is returned as the last line number in this statement, unless the entire line is a comment,
!     in which case JB is returned zero.
!  START_OF_STATEMENT is returned as the location of the first nonblank character in the statement.
!  LAST_NONBLANK is returned as the location of the last nonblank character in the statement,
!                including any trailing comments after !.
!  END_OF_STATEMENT is returned as the location of the last nonblank character in the non-comment
!                   part of the statement.

!  Any comments at the end of lines before the last line of a statement will be discarded.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: END_OF_STATEMENT,START_OF_STATEMENT,J,JB,L1,L2,L3,LAST_NONBLANK,  &
                 Q1_LEVEL,Q2_LEVEL,Q_LEVEL(132)

      JB = J
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      DO
         CALL SCAN_LINE(JB,L1,L2,L3,Q1_LEVEL,Q2_LEVEL,Q_LEVEL)
         IF (L1 == 0) THEN
             TEMP_STATEMENT(1:L3) = LINE(JB)(1:L3)
             QUOTE_LEVEL(1:L3) = Q_LEVEL(1:L3)
             START_OF_STATEMENT = L1
             END_OF_STATEMENT = L2
             LAST_NONBLANK = L3
             JB = 0
             RETURN
         ENDIF
         IF (LINE(JB)(L1:L1) == '!') THEN
             TEMP_STATEMENT(1:L3) = LINE(JB)(1:L3)
             QUOTE_LEVEL(1:L3) = Q_LEVEL(1:L3)
             START_OF_STATEMENT = L1
             END_OF_STATEMENT = L2
             LAST_NONBLANK = L3
             JB = 0
             RETURN
         ENDIF
         IF (JB == J) THEN
             TEMP_STATEMENT(1:L3) = LINE(JB)(1:L3)
             QUOTE_LEVEL(1:L3) = Q_LEVEL(1:L3)
             START_OF_STATEMENT = L1
             END_OF_STATEMENT = L2
             LAST_NONBLANK = L3
             IF (LINE(JB)(L2:L2) == '&') THEN
                 DO
                    END_OF_STATEMENT = END_OF_STATEMENT - 1
                    IF (TEMP_STATEMENT(END_OF_STATEMENT:END_OF_STATEMENT) /= ' ') EXIT
                 ENDDO
                 LAST_NONBLANK = END_OF_STATEMENT
                 JB = JB + 1
             ELSE
                 RETURN
             ENDIF
         ELSE
             TEMP_STATEMENT(LAST_NONBLANK+1:LAST_NONBLANK+L3-L1+1) = LINE(JB)(L1:L3)
             QUOTE_LEVEL(LAST_NONBLANK+1:LAST_NONBLANK+L3-L1+1) = Q_LEVEL(L1:L3)
             END_OF_STATEMENT = LAST_NONBLANK+L2-L1+1
             LAST_NONBLANK = LAST_NONBLANK+L3-L1+1
             IF (LINE(JB)(L2:L2) == '&') THEN
                 DO
                    END_OF_STATEMENT = END_OF_STATEMENT - 1
                    IF (TEMP_STATEMENT(END_OF_STATEMENT:END_OF_STATEMENT) /= ' ') EXIT
                 ENDDO
                 LAST_NONBLANK = END_OF_STATEMENT
                 JB = JB + 1
             ELSE
                 RETURN
             ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE GET_STATEMENT

      SUBROUTINE IMPLIED_DO(STATEMENT,END_OF_STATEMENT,KA,ID,JD)

!  Scan an implied do in a read list, starting at position KA -- the first comma.
!  ID is returned with the locations of the start and stop positions of the three fields:
!     loop variable,
!     starting value,
!     stopping value.
!  JD is returned with the position of the closing parenthesis of the implied do.

!  Fancier implied loops with increments other than 1 are not supported.

      IMPLICIT NONE
      CHARACTER(*) :: STATEMENT
      INTEGER :: J,JD,KA,KL,END_OF_STATEMENT,P_LEVEL,Q1_LEVEL,Q2_LEVEL,ID(6)
      LOGICAL :: IN_NAME

      P_LEVEL = 1
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      IN_NAME = .FALSE.
      ID = 0
      DO J = KA+1, END_OF_STATEMENT
         IF (STATEMENT(J:J) == ' ') CYCLE
         IF (STATEMENT(J:J) /= "=" .AND. STATEMENT(J:J) /= "," .AND. STATEMENT(J:J) /= ")") KL = J
         IF (STATEMENT(J:J) == "(" .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             P_LEVEL = P_LEVEL + 1
             CYCLE
         ENDIF
         IF (STATEMENT(J:J) == ")" .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             P_LEVEL = P_LEVEL - 1
             IF (P_LEVEL == 0) THEN
                 ID(6) = KL
                 JD = J
                 RETURN
             ENDIF
         ENDIF
         IF (STATEMENT(J:J) == "'" .AND. Q2_LEVEL == 0) THEN
             Q1_LEVEL = 1 - Q1_LEVEL
             CYCLE
         ENDIF
         IF (STATEMENT(J:J) == '"' .AND. Q1_LEVEL == 0) THEN
             Q2_LEVEL = 1 - Q2_LEVEL
             CYCLE
         ENDIF

         IF ((.NOT. IN_NAME) .AND. ID(1) == 0 .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IN_NAME = .TRUE.
             ID(1) = J
         ELSE IF (IN_NAME .AND. STATEMENT(J:J) == "=" .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IN_NAME = .FALSE.
             ID(2) = KL
         ELSE IF ((.NOT. IN_NAME) .AND. ID(2) > 0 .AND. ID(3) == 0 .AND.  &
                  Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IN_NAME = .TRUE.
             ID(3) = J
         ELSE IF (IN_NAME .AND. STATEMENT(J:J) == "," .AND. Q1_LEVEL == 0  &
                  .AND. Q2_LEVEL == 0 .AND. P_LEVEL == 1) THEN
             IN_NAME = .FALSE.
             ID(4) = KL
         ELSE IF ((.NOT. IN_NAME) .AND. ID(4) > 0 .AND. ID(5) == 0 .AND.  &
                  Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IN_NAME = .TRUE.
             ID(5) = J
         ENDIF
      ENDDO

      DO J = 1, 6
         IF (ID(J) == 0) THEN
             ID = 0
         ENDIF
      ENDDO

      END SUBROUTINE IMPLIED_DO

      FUNCTION INDEX_ULC(STRING,SUBSTRING)

!  This is the same as the intrinsic index function, except that this version
!  makes no distinction between upper and lower case letters.

!  INDEX_ULC( '      Program C2FM' , 'PROGRAM' )  will return 7.

      IMPLICIT NONE
      CHARACTER(*) :: STRING,SUBSTRING
      CHARACTER(9999), SAVE :: ST1,ST2
      INTEGER :: INDEX_ULC,J,L1,L2

!             Convert both strings to upper case and then use INDEX.

      L1 = LEN(STRING)
      L2 = LEN(SUBSTRING)
      IF (SUBSTRING(L2:L2) /= ' ') L1 = LEN_TRIM(STRING)
      DO J = 1, L1
         IF (IACHAR(STRING(J:J)) >= 97 .AND. IACHAR(STRING(J:J)) <= 122) THEN
             ST1(J:J) = ACHAR(IACHAR(STRING(J:J))-32)
         ELSE
             ST1(J:J) = STRING(J:J)
         ENDIF
      ENDDO
      DO J = 1, L2
         IF (IACHAR(SUBSTRING(J:J)) >= 97 .AND. IACHAR(SUBSTRING(J:J)) <= 122) THEN
             ST2(J:J) = ACHAR(IACHAR(SUBSTRING(J:J))-32)
         ELSE
             ST2(J:J) = SUBSTRING(J:J)
         ENDIF
      ENDDO

      INDEX_ULC = INDEX(ST1(1:L1),ST2(1:L2))

      END FUNCTION INDEX_ULC

      FUNCTION LAST_EXECUTABLE(LAST)

!  Find the last executable statement in the current routine.
!  Some new statements may need to be inserted there.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132), SAVE :: TEMP
      INTEGER :: LAST_EXECUTABLE,FIRST_NONBLANK,J,K,K1,K2,KA,LAST,INDEX_ULC
      LOGICAL :: CONTINUATION

      CONTINUATION = .FALSE.
      LAST_EXECUTABLE = 0
      DO J = 1, LAST-1
         TEMP = LINE(J)

         K1 = FIRST_NONBLANK(TEMP)
         IF (K1 < 1) CYCLE

!             Check for continuation lines.

         KA = 1
         DO K = 1, 132
            IF (TEMP(K:K) == '!') EXIT
            IF (TEMP(K:K) /= ' ') KA = K
         ENDDO
         IF (TEMP(1:1) == '!') CYCLE
         IF (TEMP(1:KA) == ' ') CYCLE
         IF (.NOT. CONTINUATION) THEN
             K2 = INDEX_ULC(TEMP(1:KA),'FORMAT')
             IF (K2 > 0) THEN
                 DO K = 1, K2
                    IF (INDEX(' 0123456789',TEMP(K:K)) > 0) CYCLE
                    LAST_EXECUTABLE = J
                 ENDDO
                 CYCLE
             ENDIF
             LAST_EXECUTABLE = J
         ENDIF
         IF (TEMP(KA:KA) == '&') THEN
             CONTINUATION = .TRUE.
         ELSE
             CONTINUATION = .FALSE.
         ENDIF
      ENDDO

      END FUNCTION LAST_EXECUTABLE

      SUBROUTINE LOCAL_VARS(KR)

!  Read on unit KR the information about local multiple precision variables for the current routine
!  from C2FM.VARS and put it into arrays ROUTINE_NAME, ROUTINE_TYPE, ROUTINE_VARIABLE,
!  VARIABLE_TYPE, ARRAY_SIZE, VARIABLE_LOCATION, and VARIABLE_INITIALIZATION.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: KR
      CHARACTER(132), SAVE :: VAR_LINE

      N_VARS = 0
      DO
         READ (KR,"(A)",END=110,ERR=110) VAR_LINE
         IF (VAR_LINE(2:8) == 'Routine') EXIT
      ENDDO

      ROUTINE_NAME = ' '
      CALL ALL_CAPS(VAR_LINE(29:60),32,ROUTINE_NAME)

      IF (VAR_LINE(18:27) == 'Program   ') ROUTINE_TYPE = 1
      IF (VAR_LINE(18:27) == 'Subroutine') ROUTINE_TYPE = 2
      IF (VAR_LINE(18:27) == 'Function  ') ROUTINE_TYPE = 3
      IF (VAR_LINE(18:27) == 'Module    ') ROUTINE_TYPE = 4

      DO
         READ (KR,"(A)",END=110) VAR_LINE
         READ (KR,"(A)",END=110) VAR_LINE
         IF (VAR_LINE(2:9) /= 'Variable') RETURN

         N_VARS = N_VARS + 1
         VARIABLE_TYPE(N_VARS) = 0
         IF (VAR_LINE(18:28) == 'Type: (FM)') VARIABLE_TYPE(N_VARS) = 1
         IF (VAR_LINE(18:28) == 'Type: (ZM)') VARIABLE_TYPE(N_VARS) = 2
         IF (VAR_LINE(18:28) == 'Type: (IM)') VARIABLE_TYPE(N_VARS) = 3

         ROUTINE_VARIABLE(N_VARS) = VAR_LINE(40:132)

         READ (KR,"(A)",END=110) VAR_LINE
         READ (VAR_LINE,"(8X,I10)") ARRAY_SIZE(N_VARS)

         VARIABLE_LOCATION(N_VARS) = VAR_LINE(21:28)

         VARIABLE_INITIALIZATION(N_VARS) = VAR_LINE(40:132)
      ENDDO

  110 RETURN

      END SUBROUTINE LOCAL_VARS

      SUBROUTINE LOGICAL_IF_RETURN(LAST)

!  Make any logical ifs of the type   IF (...) RETURN
!  into block ifs, so new code can be inserted.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132), SAVE :: TEMP
      INTEGER :: FIRST_NONBLANK,J,K,K1,K2,KA,L,LAST,INDEX_ULC
      LOGICAL :: CONTINUATION,IF_STATEMENT

      CONTINUATION = .FALSE.
      IF_STATEMENT = .FALSE.
      J = 1

  110    TEMP = LINE(J)

         K1 = FIRST_NONBLANK(TEMP)
         IF (K1 < 1) GO TO 120

         KA = 1
         DO K = 1, 132
            IF (TEMP(K:K) == '!') EXIT
            IF (TEMP(K:K) /= ' ') KA = K
         ENDDO
         IF (TEMP(1:1) == '!') GO TO 120
         IF (TEMP(1:KA) == ' ') GO TO 120
         IF (.NOT. CONTINUATION) THEN

!             Check for the start of an IF statement.

             IF (INDEX_ULC(TEMP(K1:K1+1),'IF') == 1) THEN
                 DO K = K1+2, KA
                    IF (TEMP(K:K) == ' ') CYCLE
                    IF (TEMP(K:K) == '(') THEN
                        IF_STATEMENT = .TRUE.
                        EXIT
                    ELSE
                        IF_STATEMENT = .FALSE.
                        EXIT
                    ENDIF
                 ENDDO
             ELSE
                 IF_STATEMENT = .FALSE.
             ENDIF
         ENDIF

!             Check for continuation lines.

         IF (TEMP(KA:KA) == '&') THEN
             CONTINUATION = .TRUE.
         ELSE
             CONTINUATION = .FALSE.
         ENDIF

!             Check for the end of an IF statement.

         IF (IF_STATEMENT) THEN
             IF (.NOT. CONTINUATION) THEN
                 K2 = KA
                 DO K = KA, K1, -1
                    IF (TEMP(K:K) /= ')') CYCLE
                    K2 = K + 1
                    EXIT
                 ENDDO
                 DO K = K2, KA
                    IF (TEMP(K:K) /= ' ') THEN
                        IF (INDEX_ULC(TEMP(K:KA),'RETURN') == 1 .AND. KA == K+5) THEN

!                           Change the logical if to a block if.

                            LAST = LAST + 2
                            DO L = LAST, J+3, -1
                               LINE(L) = LINE(L-2)
                            ENDDO
                            LINE(J)(K:KA) = 'THEN  '
                            LINE(J+1) = ' '
                            LINE(J+1)(K1:K1+5) = 'RETURN'
                            LINE(J+2) = ' '
                            LINE(J+2)(K1:K1+4) = 'ENDIF'
                        ENDIF
                        EXIT
                    ENDIF
                 ENDDO
             ENDIF
             IF_STATEMENT = .FALSE.
         ENDIF

  120 J = J + 1
      IF (J < LAST) GO TO 110

      END SUBROUTINE LOGICAL_IF_RETURN

      SUBROUTINE READ_ROUTINE(KINP,FIRST,LAST,FILE_ENDED)

!  Read the next routine on unit KINP.
!  LINE(1), ..., LINE(LAST) contains the routine.
!  LINE(FIRST) is the start of the PROGRAM, SUBROUTINE, FUNCTION, or MODULE subprogram.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132) :: TEMP
      INTEGER :: FIRST,J,K,K1,KINP,KL(4),L,LAST,INDEX_ULC
      LOGICAL :: FILE_ENDED

      FILE_ENDED = .FALSE.
      ROUTINE_TYPE = 0
      FIRST = 0
      DO J = 1, 10000
         READ (KINP,"(A)",END=110) LINE(J)

!             Look for the first statement.

         KL(1) = INDEX_ULC(LINE(J),'PROGRAM')
         KL(2) = INDEX_ULC(LINE(J),'SUBROUTINE')
         KL(3) = INDEX_ULC(LINE(J),'FUNCTION')
         KL(4) = INDEX_ULC(LINE(J),'MODULE')
         K = 0
         DO K1 = 1, 4
            IF (K == 0 .AND. KL(K1) > 0) K = KL(K1)
            IF (K > 0  .AND. KL(K1) > 0) K = MIN(K,KL(K1))
         ENDDO

         IF (FIRST == 0) THEN
             IF (KL(1) > 0 .AND. KL(1) == K) CALL FIRST_STATEMENT(LINE(J),FIRST,J,KL(1))
             IF (FIRST > 0) ROUTINE_TYPE = 1
         ENDIF
         IF (FIRST == 0) THEN
             IF (KL(2) > 0 .AND. KL(2) == K) CALL FIRST_STATEMENT(LINE(J),FIRST,J,KL(2))
             IF (FIRST > 0) ROUTINE_TYPE = 2
         ENDIF
         IF (FIRST == 0) THEN
             IF (KL(3) > 0 .AND. KL(3) == K) CALL FIRST_STATEMENT(LINE(J),FIRST,J,KL(3))
             IF (FIRST > 0) ROUTINE_TYPE = 3
         ENDIF
         IF (FIRST == 0) THEN
             IF (KL(4) > 0 .AND. KL(4) == K) CALL FIRST_STATEMENT(LINE(J),FIRST,J,KL(4))
             IF (FIRST > 0) ROUTINE_TYPE = 4
         ENDIF

!             Look for the last statement.

         IF (FIRST == 0) CYCLE
         K = INDEX_ULC(LINE(J),'END ')
         IF (K > 0) THEN
             TEMP = ' '
             IF (K > 1) THEN
                 IF (LINE(J)(1:K-1) /= TEMP(1:K-1)) CYCLE
             ENDIF
             IF (K+3 > 132) THEN
                 LAST = J
                 RETURN
             ENDIF
             IF (LINE(J)(K+3:132) == TEMP(K+3:132)) THEN
                 LAST = J
                 RETURN
             ENDIF
             L = INDEX(LINE(J),'!')
             K = INDEX_ULC(LINE(J),'PROGRAM')
             IF (K > 0 .AND. (L == 0 .OR. L > K)) THEN
                 LAST = J
                 RETURN
             ENDIF
             K = INDEX_ULC(LINE(J),'SUBROUTINE')
             IF (K > 0 .AND. (L == 0 .OR. L > K)) THEN
                 LAST = J
                 RETURN
             ENDIF
             K = INDEX_ULC(LINE(J),'FUNCTION')
             IF (K > 0 .AND. (L == 0 .OR. L > K)) THEN
                 LAST = J
                 RETURN
             ENDIF
             K = INDEX_ULC(LINE(J),'MODULE')
             IF (K > 0 .AND. (L == 0 .OR. L > K)) THEN
                 LAST = J
                 RETURN
             ENDIF
         ENDIF

!             Internal subprograms will be scanned the same way as external ones.
!             Treat a CONTAINS statement as the end of the current routine.

         K = INDEX_ULC(LINE(J),'CONTAINS')
         IF (K > 0) THEN
             TEMP = ' '
             IF (K > 1) THEN
                 IF (LINE(J)(1:K-1) /= TEMP(1:K-1)) CYCLE
             ENDIF
             IF (K+8 > 132) THEN
                 LAST = J
                 RETURN
             ENDIF
             IF (LINE(J)(K+8:132) == TEMP(K+8:132)) THEN
                 LAST = J
                 RETURN
             ENDIF
         ENDIF

         CYCLE

!             End of file encountered.

  110    LAST = J - 1
         FILE_ENDED = .TRUE.
         RETURN
      ENDDO

      END SUBROUTINE READ_ROUTINE

      SUBROUTINE READ_VAR_SUBSCRIPT(READ_VAR)

!  Fix the constants in newly generated variables in read statements.
!  For example, READ_VAR might be C2FM_TEMP(1+L+1+1).
!  This routine converts that to a more readable form: C2FM_TEMP(L+3).

      IMPLICIT NONE
      CHARACTER(64) :: READ_VAR,TEMP
      CHARACTER(10) :: NUMB
      INTEGER :: ITEM,J,JA,JB,JA_START,JT,K,Q1_LEVEL,Q2_LEVEL,P_LEVEL,VALUE,SUM_VALUE
      LOGICAL :: NON_CONSTANT

      TEMP = ' '
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      P_LEVEL = 0
      JT = 0
      NON_CONSTANT = .FALSE.
      ITEM = 0
      VALUE = 0
      SUM_VALUE = 0

      JA = 1
      JB = 0
      JA_START = 1

  110 DO J = JA+1, 64
         IF (READ_VAR(J:J) == ' ') CYCLE
         IF (READ_VAR(J:J) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) P_LEVEL = P_LEVEL + 1
         IF (READ_VAR(J:J) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (READ_VAR(J:J) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (Q1_LEVEL > 0 .OR. Q2_LEVEL > 0 .OR. P_LEVEL /= 1) CYCLE
         IF (READ_VAR(J:J) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) P_LEVEL = P_LEVEL - 1

         IF (INDEX('(+-:,)',READ_VAR(J:J)) > 0) THEN
             IF (JA == JA_START) THEN
                 JA = J
                 TEMP(JT+1:JT+JA-JA_START+1) = READ_VAR(JA_START:JA)
                 JT = JT + JA - JA_START + 1
                 GO TO 110
             ELSE
                 JB = J
                 GO TO 120
             ENDIF
         ENDIF

         IF (INDEX('0123456789',READ_VAR(J:J)) == 0) NON_CONSTANT = .TRUE.
         CYCLE

!             If the item between the last two delimiters is a constant, convert it.

  120    IF (NON_CONSTANT) THEN
             ITEM = ITEM + 1
             IF (ITEM == 1) THEN
                 TEMP(JT+1:JT+JB-JA-1) = READ_VAR(JA+1:JB-1)
                 JT = JT + JB - JA - 1
                 JA = JB
             ELSE
                 TEMP(JT+1:JT+JB-JA) = READ_VAR(JA:JB-1)
                 JT = JT + JB - JA
                 JA = JB
             ENDIF
         ELSE
             READ (READ_VAR(JA+1:JB-1),*) VALUE
             IF (READ_VAR(JA:JA) == '-') THEN
                 SUM_VALUE = SUM_VALUE - VALUE
             ELSE
                 SUM_VALUE = SUM_VALUE + VALUE
             ENDIF
         ENDIF
         NON_CONSTANT = .FALSE.

!             See if this is the end of a subscript.

         IF (INDEX(':,)',READ_VAR(J:J)) > 0) THEN
             IF (SUM_VALUE > 0) THEN
                 WRITE (NUMB,"(I10)") SUM_VALUE
                 IF (ITEM > 0) THEN
                     JT = JT + 1
                     TEMP(JT:JT) = '+'
                 ENDIF
                 DO K = 1, 10
                    IF (NUMB(K:K) /= ' ') THEN
                        JT = JT + 1
                        TEMP(JT:JT) = NUMB(K:K)
                    ENDIF
                 ENDDO
             ELSE IF (SUM_VALUE < 0) THEN
                 WRITE (NUMB,"(I10)") SUM_VALUE
                 DO K = 1, 10
                    IF (NUMB(K:K) /= ' ') THEN
                        JT = JT + 1
                        TEMP(JT:JT) = NUMB(K:K)
                    ENDIF
                 ENDDO
             ENDIF
             IF (READ_VAR(J:J) == ')') THEN
                 JT = JT + 1
                 TEMP(JT:JT) = ')'
                 EXIT
             ENDIF
             JT = JT + 1
             TEMP(JT:JT) = READ_VAR(J:J)
             SUM_VALUE = 0
             ITEM = 0
             JA = JB
             JB = JA
             JA_START = JA - 1
             NON_CONSTANT = .FALSE.
             GO TO 110
         ELSE
             JA = JB
             GO TO 110
         ENDIF
      ENDDO

      READ_VAR = TEMP

      END SUBROUTINE READ_VAR_SUBSCRIPT

      SUBROUTINE SAVE_LOCAL_VARS

!  If there are any multiple precision variables that are local to a subroutine or function,
!  Add a SAVE statement for them.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: J,K,L

      IF (ROUTINE_TYPE == 1) RETURN

      TEMP_STATEMENT = '      SAVE :: '

      K = 15
      DO J = 1, N_VARS
         IF (VARIABLE_TYPE(J) == 0) CYCLE
         IF (VARIABLE_LOCATION(J) == 'Local   ') THEN
             L = LEN_TRIM(ROUTINE_VARIABLE(J))
             IF (K == 15) THEN
                 TEMP_STATEMENT(K:K+L-1) = ROUTINE_VARIABLE(J)(1:L)
                 K = K + L
             ELSE
                 TEMP_STATEMENT(K:K+L+1) = ', '//ROUTINE_VARIABLE(J)(1:L)
                 K = K + L + 2
             ENDIF
         ENDIF
      ENDDO

      K = K - 1
      J = K
      IF (K > 14) CALL WRITE_STATEMENT(J,K)

      END SUBROUTINE SAVE_LOCAL_VARS

      SUBROUTINE SCAN_LINE(J,L1,L2,L3,Q1_LEVEL,Q2_LEVEL,Q_LEVEL)

!  LINE(J) contains part of the next statement.
!  L1 is returned as the first nonblank character on the line.
!  L2 is returned as the last nonblank character on the line, except for comments.
!  L3 is returned as the last nonblank character on the line.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: J,K,L1,L2,L3,COMMENT,Q1_LEVEL,Q2_LEVEL,Q_LEVEL(132)

      COMMENT = 0
      L1 = 0
      L2 = 0
      L3 = 0
      DO K = 1, 132
         IF (L1 == 0 .AND. LINE(J)(K:K) /= ' ') L1 = K
         IF (LINE(J)(K:K) == "'" .AND. Q2_LEVEL == 0 .AND. COMMENT == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (LINE(J)(K:K) == '"' .AND. Q1_LEVEL == 0 .AND. COMMENT == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND. LINE(J)(K:K) == '!') COMMENT = 1
         Q_LEVEL(K) = MAX(Q1_LEVEL,Q2_LEVEL)
         IF (COMMENT == 0 .AND. LINE(J)(K:K) /= ' ') L2 = K
         IF (LINE(J)(K:K) /= ' ') L3 = K
      ENDDO

      END SUBROUTINE SCAN_LINE

      SUBROUTINE SCAN_LIST(STATEMENT,KA,END_OF_STATEMENT,KB)

!  STATEMENT contains the statement.
!  Scan from position KA to END_OF_STATEMENT looking for an argument list or
!  subscript expression.
!
!  KB is returned with the position of the final ) in the list, or zero if there is
!     no list.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(*) :: STATEMENT
      INTEGER :: J,KA,KB,END_OF_STATEMENT,P_LEVEL,Q1_LEVEL,Q2_LEVEL
      LOGICAL :: IN_LIST

      KB = 0
      P_LEVEL = 0
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      IN_LIST = .FALSE.
      DO J = KA, END_OF_STATEMENT
         IF (STATEMENT(J:J) == ' ') CYCLE
         IF ((.NOT. IN_LIST) .AND. STATEMENT(J:J) /= '(') CYCLE
         IF ((.NOT. IN_LIST) .AND. STATEMENT(J:J) == '(') IN_LIST = .TRUE.
         IF (STATEMENT(J:J) == "(" .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             P_LEVEL = P_LEVEL + 1
         ENDIF
         IF (STATEMENT(J:J) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (STATEMENT(J:J) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (STATEMENT(J:J) == ")" .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             P_LEVEL = P_LEVEL - 1
             IF (P_LEVEL == 0) THEN
                 KB = J
                 RETURN
             ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE SCAN_LIST

      SUBROUTINE SET_TYPE(KTYPE,L,TOTAL_VARS)

!  LINE(L) is the start of a type declaration statement in which the variables will be
!          converted to multiple precision.

!  KTYPE = 1, 2, 3 for conversion to types (FM), (ZM), (IM).

!  ROUTINE_ARGUMENT(1:N_ARGS) are the input arguments to this routine, and now the type
!                             declarations are scanned to identify in ARGUMENT_TYPE which
!                             ones are being converted to multiple precision.

!  ARGUMENT_TYPE(K) = 1 for TYPE (FM)
!                   = 2 for TYPE (ZM)
!                   = 3 for TYPE (IM)

!  Also, local multiple precision variables that are not input arguments are identified.
!  Information is written to C2FM.VARS about all multiple precision variables used in each routine.

!  N_VARS is the number of variables in the current routine.

!  ROUTINE_VARIABLE(K) is the variable name of the Kth variable in the current routine that
!                      will become multiple precision.

!  VARIABLE_INITIALIZATION(K) gives any initialization expression in the declaration of
!                             the original program.  These must be made into executable
!                             statements and put at the top of the routine.

!  VARIABLE_LOCATION(K) records whether that variable is an input argument, a module variable,
!                       or a local variable.

!  ARRAY_SIZE(K) gives the size of a multiple precision variable array.

!  VARIABLE_TYPE(K) = 0, 1, 2, 3 gives these local variable types as with ARGUMENT_TYPE.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(32)  :: TEMP
      INTEGER :: FIRST,K,KA,KB,K1,K2,KTYPE,L,L2,TOTAL_VARS,  &
                 Q1_LEVEL,Q2_LEVEL,P_LEVEL,DEFAULT_SIZE,GET_SIZE,INDEX_ULC

!             Scan the variables in the declaration.

!             Q1_LEVEL or Q2_LEVEL is positive when the scan in within a character string.
!             P_LEVEL is positive when the scan in within a pair of parentheses.

      FIRST = L
      Q1_LEVEL = 0
      Q2_LEVEL = 0
      P_LEVEL = 0
      N_VARS = 0
      ARRAY_SIZE = 1
      DEFAULT_SIZE = 1
      VARIABLE_LOCATION = 'Local   '
      VARIABLE_INITIALIZATION = ' '
      ROUTINE_VARIABLE(1) = ' '
      VARIABLE_TYPE(1:100) = -1

!             Start the variable scan at position KA.

!             If there is a DIMENSION attribute present in the declaration, that becomes
!             the default size.

      KA = INDEX(LINE(FIRST),'::')
      IF (KA <= 0) THEN
          IF (KTYPE <= 2) THEN
              KA = INDEX_ULC(LINE(FIRST),'(KIND(3.1D1))')
              KA = KA + 13
          ELSE
              KA = INDEX_ULC(LINE(FIRST),'(KIND(31))')
              KA = KA + 10
          ENDIF
      ELSE
          IF (INDEX_ULC(LINE(FIRST)(1:KA),'ALLOCATABLE') > 0) THEN
              ARRAY_SIZE = 0
              DEFAULT_SIZE = 0
          ENDIF
          IF (INDEX_ULC(LINE(FIRST)(1:KA),'EXTERNAL') > 0) VARIABLE_LOCATION = 'External'
          KB = INDEX_ULC(LINE(FIRST)(1:KA),'DIMENSION')
          IF (KB <= 0) THEN
              KA = KA + 2
              GO TO 130
          ENDIF
          K1 = KA
          DO L2 = KB, K1
             KA = L2
             IF (LINE(FIRST)(KA:KA) == '(') EXIT
          ENDDO
          IF (LINE(FIRST)(KA:KA) == '(') THEN
              P_LEVEL = 1
  110         KB = KA + 1
              DO L2 = KB, 132
                 IF (LINE(FIRST)(L2:L2) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                     P_LEVEL = P_LEVEL - 1
                     IF (P_LEVEL == 0) THEN
                         IF (DEFAULT_SIZE > 0) THEN
                             DEFAULT_SIZE = GET_SIZE(LINE(FIRST),KA,L2)
                         ENDIF
                         KA = L2 + 1
                         EXIT
                     ENDIF
                 ENDIF
                 IF (LINE(FIRST)(L2:L2) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                 IF (LINE(FIRST)(L2:L2) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                 IF (LINE(FIRST)(L2:L2) == '&' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
  120                FIRST = FIRST + 1
                     DO K = 1, 132
                        IF (LINE(FIRST)(K:K) /= ' ') THEN
                            IF (LINE(FIRST)(K:K) == '!') THEN
                                GO TO 120
                            ELSE
                                KA = K - 1
                                GO TO 110
                            ENDIF
                        ENDIF
                     ENDDO
                     GO TO 120
                 ENDIF
              ENDDO
          ENDIF
          ARRAY_SIZE = DEFAULT_SIZE
          KA = K1 + 2
      ENDIF

!             Advance KA to the next non-blank character, to find the start of the next name.

  130 K1 = 0
      DO K = KA, 132
         IF (LINE(FIRST)(K:K) /= ' ') THEN
             IF (LINE(FIRST)(K:K) /= '!') K1 = K
             EXIT
         ENDIF
      ENDDO
      IF (K1 == 0) GO TO 200
      KA = K1
      IF (LINE(FIRST)(K:K) == '&') THEN
          FIRST = FIRST + 1
          KA = 1
          GO TO 130
      ENDIF

!             Scan the name.

      KB = K1
      DO L2 = KB, 132
         IF (LINE(FIRST)(L2:L2) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             P_LEVEL = P_LEVEL + 1
             IF (P_LEVEL == 1) THEN
                 KA = L2
                 N_VARS = N_VARS + 1
                 ROUTINE_VARIABLE(N_VARS) = LINE(FIRST)(K1:L2-1)
                 GO TO 150
             ENDIF
         ENDIF
         IF (LINE(FIRST)(L2:L2) == ' ' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IF (P_LEVEL == 0) THEN
                 KA = L2
                 N_VARS = N_VARS + 1
                 ROUTINE_VARIABLE(N_VARS) = LINE(FIRST)(K1:L2-1)
                 GO TO 150
             ENDIF
         ENDIF
         IF (LINE(FIRST)(L2:L2) == '=' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IF (P_LEVEL == 0) THEN
                 KA = L2
                 N_VARS = N_VARS + 1
                 ROUTINE_VARIABLE(N_VARS) = LINE(FIRST)(K1:L2-1)
                 GO TO 180
             ENDIF
         ENDIF
         IF (LINE(FIRST)(L2:L2) == ',' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
             IF (P_LEVEL == 0) THEN
                 KA = L2 + 1
                 N_VARS = N_VARS + 1
                 ROUTINE_VARIABLE(N_VARS) = LINE(FIRST)(K1:L2-1)
                 GO TO 130
             ENDIF
         ENDIF
         IF (LINE(FIRST)(L2:L2) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
         IF (LINE(FIRST)(L2:L2) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
         IF (LINE(FIRST)(L2:L2) == '&' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
             P_LEVEL == 0) THEN
  140        FIRST = FIRST + 1
             DO K = 1, 132
                IF (LINE(FIRST)(K:K) /= ' ') THEN
                    IF (LINE(FIRST)(K:K) == '!') THEN
                        GO TO 140
                    ELSE
                        KA = K
                        GO TO 130
                    ENDIF
                ENDIF
             ENDDO
             GO TO 140
         ENDIF
      ENDDO
      GO TO 200

!             Advance KA to the next non-blank character.

  150 K1 = 0
      DO K = KA, 132
         IF (LINE(FIRST)(K:K) /= ' ') THEN
             K1 = K
             EXIT
         ENDIF
      ENDDO
      IF (K1 == 0) GO TO 200
      KA = K1

!             Scan any dimension information.

      IF (LINE(FIRST)(KA:KA) == '(') THEN
  160     KB = KA + 1
          DO L2 = KB, 132
             IF (LINE(FIRST)(L2:L2) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
                 P_LEVEL = P_LEVEL - 1
                 IF (P_LEVEL == 0) THEN
                     IF (DEFAULT_SIZE > 0) THEN
                         ARRAY_SIZE(N_VARS) = GET_SIZE(LINE(FIRST),KA,L2)
                     ENDIF
                     KA = L2 + 1
                     EXIT
                 ENDIF
             ENDIF
             IF (LINE(FIRST)(L2:L2) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
             IF (LINE(FIRST)(L2:L2) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
             IF (LINE(FIRST)(L2:L2) == '&' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0) THEN
  170            FIRST = FIRST + 1
                 DO K = 1, 132
                    IF (LINE(FIRST)(K:K) /= ' ') THEN
                        IF (LINE(FIRST)(K:K) == '!') THEN
                            GO TO 170
                        ELSE
                            KA = K - 1
                            GO TO 160
                        ENDIF
                    ENDIF
                 ENDDO
                 GO TO 170
             ENDIF
          ENDDO
      ENDIF

!             Scan any initialization expression.  (Must be entirely on this line)

  180 K1 = KA
      DO L2 = K1, 132
         IF (LINE(FIRST)(L2:L2) /= ' ') THEN
             IF (LINE(FIRST)(L2:L2) == ',') THEN
                 KA = KA + 1
                 GO TO 130
             ENDIF
             IF (LINE(FIRST)(L2:L2) == '=') THEN
                 K2 = L2
                 KA = L2 + 1
                 DO K = L2, 132
                    IF (LINE(FIRST)(K:K) == '"' .AND. Q1_LEVEL == 0) Q2_LEVEL = 1 - Q2_LEVEL
                    IF (LINE(FIRST)(K:K) == "'" .AND. Q2_LEVEL == 0) Q1_LEVEL = 1 - Q1_LEVEL
                    IF (LINE(FIRST)(K:K) == '(' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0)  &
                        P_LEVEL = P_LEVEL + 1
                    IF (LINE(FIRST)(K:K) == ')' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0)  &
                        P_LEVEL = P_LEVEL - 1
                    IF (LINE(FIRST)(K:K) == ',' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
                        P_LEVEL == 0) THEN
                        KB = K - 1
                        GO TO 190
                    ENDIF
                    IF (LINE(FIRST)(K:K) == '&' .AND. Q1_LEVEL == 0 .AND. Q2_LEVEL == 0 .AND.  &
                        P_LEVEL == 0) THEN
                        KB = K - 1
                        GO TO 190
                    ENDIF
                 ENDDO
                 KB = 132
                 GO TO 190
             ENDIF
         ENDIF
      ENDDO
      GO TO 130

  190 VARIABLE_INITIALIZATION(N_VARS) = TRIM(LINE(FIRST)(KA:KB))
      DO K = KA, KB
         IF (LINE(FIRST)(K:K) /= ' ') THEN
             VARIABLE_INITIALIZATION(N_VARS) = TRIM(LINE(FIRST)(K:KB))
             EXIT
         ENDIF
      ENDDO
      LINE(FIRST)(K2:KB) = ' '
      IF (KB < 131) THEN
          KA = KB + 2
          IF (LINE(FIRST)(KB+1:KB+1) == '&') KA = KB + 1
          GO TO 130
      ENDIF

!             The end of this declaration statement has been reached.

!             Convert all names to upper case to make later comparisons easier.

  200 DO K = 1, N_ARGS
         CALL ALL_CAPS(ROUTINE_ARGUMENT(K),32,TEMP)
         ROUTINE_ARGUMENT(K) = TEMP
      ENDDO
      DO K = 1, N_VARS
         CALL ALL_CAPS(ROUTINE_VARIABLE(K),32,TEMP)
         ROUTINE_VARIABLE(K) = TEMP
      ENDDO

!             Check if this multiple precision declaration includes any variables in this
!             routine's argument list.  If so, set ARGUMENT_TYPE.

      DO K = 1, N_ARGS
         DO L2 = 1, N_VARS
            IF (ROUTINE_ARGUMENT(K) == ROUTINE_VARIABLE(L2)) THEN
                ARGUMENT_TYPE(K) = KTYPE
                VARIABLE_LOCATION(L2) = 'Argument'
            ENDIF
         ENDDO
      ENDDO

!             Check if this multiple precision declaration includes the current function
!             name.  If so, set VARIABLE_LOCATION.

      DO L2 = 1, N_VARS
         IF (ROUTINE_TYPE == 3 .AND. ROUTINE_NAME == ROUTINE_VARIABLE(L2)) THEN
             VARIABLE_LOCATION(L2) = 'External'
         ENDIF
      ENDDO

!             Write the variable information to file C2FM.VARS for ordinary routines, also
!             to file C2FM.MODS for modules.

      IF (ROUTINE_TYPE == 4) THEN
          L2 = 33
      ELSE
          L2 = 32
      ENDIF

      DO KA = 32, L2
         DO K = 1, N_VARS
            TOTAL_VARS = TOTAL_VARS + 1
            IF (KTYPE == 1) THEN
                WRITE (KA,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (FM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (KA,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ELSE IF (KTYPE == 2) THEN
                WRITE (KA,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (ZM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (KA,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ELSE IF (KTYPE == 3) THEN
                WRITE (KA,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (IM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (KA,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE SET_TYPE

      SUBROUTINE SORT_NAMES

!  Sort the R_NAME database by routine names to speed up later searches.

!  Shell sort.
!  The increments used are:  NH, ..., 40, 13, 4, 1, where
!  NH is defined by 3*NH + 1 < N_ROUTINES < 3*(3*NH + 1) + 1

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(32) :: RT
      INTEGER :: I,IN,J,JA,K,NH,NH1,NH2,RU,RV,RW,RX,RY

!            Convert the names to upper case.

      DO J = 1, N_ROUTINES
         CALL ALL_CAPS(R_NAME(J),32,RT)
         R_NAME(J) = RT
      ENDDO

!            Find NH for pass 1.

      IF (N_ROUTINES <= 1) RETURN
      NH2 = 4
      NH1 = 1
  110 NH = NH1
      NH1 = NH2
      NH2 = 3*NH2 + 1
      IF (NH2 < N_ROUTINES) GO TO 110

!            Begin the next pass.

  120 DO J = 1, NH

!            Insertion sort on the list R_NAME(J), R_NAME(J+NH), ..., R_NAME(J+K*NH)

         K = (N_ROUTINES-J)/NH
         DO JA = 1, K

!            Insert R_NAME(J+JA*NH) into the sorted list
!            R_NAME(J), ..., R_NAME(J+(JA-1)*NH).

            IN = J + JA*NH
            I = IN - NH
            RT = R_NAME(IN)
            RU = R_TYPE(IN)
            RV = R_N_ARGS(IN)
            RW = R_ARG_START(IN)
            RX = R_N_USED(IN)
            RY = R_USED_START(IN)

  130       IF (RT >= R_NAME(I)) GO TO 140
            R_NAME(IN) = R_NAME(I)
            R_TYPE(IN) = R_TYPE(I)
            R_N_ARGS(IN) = R_N_ARGS(I)
            R_ARG_START(IN) = R_ARG_START(I)
            R_N_USED(IN) = R_N_USED(I)
            R_USED_START(IN) = R_USED_START(I)

            IN = I
            I = I - NH
            IF (I > 0) GO TO 130

  140       R_NAME(IN) = RT
            R_TYPE(IN) = RU
            R_N_ARGS(IN) = RV
            R_ARG_START(IN) = RW
            R_N_USED(IN) = RX
            R_USED_START(IN) = RY
         ENDDO
      ENDDO

      NH = NH/3
      IF (NH > 0) GO TO 120

      END SUBROUTINE SORT_NAMES

      SUBROUTINE USE_CHECK

!  Expand the R_NAME database by adding any modules used by other modules to
!  the list of used modules in each routine.
!  This creates the R_N_ALL_USED and R_ALL_USED part of the database.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: J,J2,K,N,N_MODS,NEW_TOTAL,TOTAL_USED

!             Count the number of modules present in the program.

      N_MODS = 0
      DO J = 1, N_ROUTINES
         IF (R_TYPE(J) == 4) N_MODS = N_MODS + 1
         R_N_ALL_USED(J) = R_N_USED(J)
      ENDDO
      ALLOCATE ( R_ALL_USED(N_ROUTINES,N_MODS) )
      R_ALL_USED = 0

!             Initialize R_N_ALL_USED with the modules used directly.

      DO J = 1, N_ROUTINES
         IF (R_N_USED(J) > 0) THEN
             DO K = 1, R_N_USED(J)
                CALL FIND_NAME(R_USED(R_USED_START(J)+K-1),R_NAME,N_ROUTINES,J2)
                R_ALL_USED(J,K) = J2
             ENDDO
         ENDIF
      ENDDO

      DO
         TOTAL_USED = 0
         DO J = 1, N_ROUTINES
            TOTAL_USED = TOTAL_USED + R_N_ALL_USED(J)
         ENDDO
         DO J = 1, N_ROUTINES
            IF (R_N_ALL_USED(J) > 0) THEN
                N = R_N_ALL_USED(J)
                DO K = 1, N
                   CALL EXPAND_USE_LIST(J,R_ALL_USED(J,K))
                ENDDO
            ENDIF
         ENDDO
         NEW_TOTAL = 0
         DO J = 1, N_ROUTINES
            NEW_TOTAL = NEW_TOTAL + R_N_ALL_USED(J)
         ENDDO
         IF (NEW_TOTAL == TOTAL_USED) EXIT
      ENDDO

      END SUBROUTINE USE_CHECK

      SUBROUTINE VAR2_LIST

!  Make a second variable list for each routine, including any module variables that are used.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(132), SAVE :: VAR_LINE
      INTEGER :: I,J,K,L,TOTAL_VARS

      REWIND(32)
      DO J = 1, N_ROUTINES

!             Get the information about local multiple precision variables in routine J.

         CALL LOCAL_VARS(32)

!             Routine J in the file is number K in the sorted database.  Find K

         K = 0
         DO L = 1, N_ROUTINES
            IF (ROUTINE_NAME == R_NAME(L)) THEN
                K = L
                EXIT
            ENDIF
         ENDDO

!             Get the module variables for each used module.

         IF (R_N_ALL_USED(K) > 0) THEN
             DO I = 1, R_N_ALL_USED(K)
                REWIND(33)
                DO
                   READ (33,"(A)",END=110) VAR_LINE
                   IF (VAR_LINE(2:8) == 'Routine' .AND.  &
                       VAR_LINE(29:132) == R_NAME(R_ALL_USED(K,I))) EXIT
                ENDDO
                DO
                   READ (33,"(A)",END=110) VAR_LINE
                   READ (33,"(A)",END=110) VAR_LINE
                   IF (VAR_LINE(2:9) /= 'Variable') EXIT
                   N_VARS = N_VARS + 1
                   VARIABLE_TYPE(N_VARS) = 0
                   IF (VAR_LINE(18:28) == 'Type: (FM)') VARIABLE_TYPE(N_VARS) = 1
                   IF (VAR_LINE(18:28) == 'Type: (ZM)') VARIABLE_TYPE(N_VARS) = 2
                   IF (VAR_LINE(18:28) == 'Type: (IM)') VARIABLE_TYPE(N_VARS) = 3
                   ROUTINE_VARIABLE(N_VARS) = VAR_LINE(40:132)
                   READ (33,"(A)",END=110) VAR_LINE
                   READ (VAR_LINE,"(8X,I10)") ARRAY_SIZE(N_VARS)
                   VARIABLE_LOCATION(N_VARS) = 'Module  '
                   VARIABLE_INITIALIZATION(N_VARS) = VAR_LINE(40:132)
                ENDDO
             ENDDO
         ENDIF

!             Write the second variable list.

  110    IF (ROUTINE_TYPE == 4) THEN
             WRITE (35,"(//A,I4,A,A/)") ' Routine ',J,':   Module     ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 3) THEN
             WRITE (35,"(//A,I4,A,A/)") ' Routine ',J,':   Function   ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 2) THEN
             WRITE (35,"(//A,I4,A,A/)") ' Routine ',J,':   Subroutine ',TRIM(ROUTINE_NAME)
         ELSE IF (ROUTINE_TYPE == 1) THEN
             WRITE (35,"(//A,I4,A,A/)") ' Routine ',J,':   Program    ',TRIM(ROUTINE_NAME)
         ENDIF
         TOTAL_VARS = 0
         DO K = 1, N_VARS
            TOTAL_VARS = TOTAL_VARS + 1
            IF (VARIABLE_TYPE(K) == 1) THEN
                WRITE (35,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (FM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (35,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ELSE IF (VARIABLE_TYPE(K) == 2) THEN
                WRITE (35,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (ZM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (35,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ELSE IF (VARIABLE_TYPE(K) == 3) THEN
                WRITE (35,"(A,I4,A,A)") ' Variable ',TOTAL_VARS,'   Type: (IM)  ,  Name = ',  &
                                         TRIM(ROUTINE_VARIABLE(K))
                WRITE (35,"(A,I10,A,A,A,A/)") ' Size = ',ARRAY_SIZE(K),', ',VARIABLE_LOCATION(K),  &
                                          ',  Init =  ',TRIM(VARIABLE_INITIALIZATION(K))
            ENDIF
         ENDDO

      ENDDO

      END SUBROUTINE VAR2_LIST

      SUBROUTINE WRITE_PASS2(LAST)

!  Write the next routine.
!  LINE(1), ..., LINE(LAST) contains the routine.

      USE C2FM_VARS
      IMPLICIT NONE
      CHARACTER(32), SAVE :: V_NAME
      CHARACTER(132), SAVE :: VAR_LINE,V_INIT
      INTEGER :: FIRST_NONBLANK,J,J1,J2,JB,K,LAST,LAST_EXECUTABLE,INDEX_ULC

!             Make any logical ifs of the type   IF (...) RETURN
!             into block ifs, so new code can be inserted.

      CALL LOGICAL_IF_RETURN(LAST)

      CALL FIRST_EXECUTABLE(LAST,J1,J2)

      DO J = 1, LAST
         IF (J == J1 .AND. N_VARS > 0) THEN
             WRITE (36,"(A)") '      USE FMZM'
             WRITE (36,"(A)") '      USE C2FM_READS'
         ENDIF

         IF (J == J2) THEN

!             See if any multiple precision variables were originally initialized.

             JB = 0
             DO K = 1, N_VARS
                IF (VARIABLE_INITIALIZATION(K) /= ' ' .AND. ROUTINE_TYPE /= 4) THEN
                    JB = JB + 1
                    IF (JB == 1) THEN
                        WRITE (36,"(A)") ' '
                        IF (ROUTINE_TYPE /= 1) THEN
                            WRITE (36,"(A)") '      LOGICAL, SAVE :: FIRST_CALL_FM = .TRUE.'
                            WRITE (36,"(A)") ' '
                        ENDIF
                        IF (ROUTINE_TYPE == 1) THEN
                            WRITE (36,"(A)") '      CALL FM_SET(50)'
                        ELSE IF (ROUTINE_TYPE == 2) THEN
                            WRITE (36,"(A)") '      CALL FM_ENTER_USER_ROUTINE'
                        ELSE IF (ROUTINE_TYPE == 3) THEN
                            WRITE (36,"(A,A,A)") '      CALL FM_ENTER_USER_FUNCTION('  &
                                             ,TRIM(ROUTINE_NAME),')'
                        ENDIF
                        WRITE (36,"(A)") ' '
                        WRITE (36,"(A,A)") '!             Multiple precision variables that ',  &
                                           'were originally initialized in declarations.'
                        WRITE (36,"(A)") ' '
                        IF (ROUTINE_TYPE /= 1) THEN
                            WRITE (36,"(A)") '      IF (FIRST_CALL_FM) THEN'
                            WRITE (36,"(A)") '          FIRST_CALL_FM = .FALSE.'
                        ENDIF
                    ENDIF
                    IF (VARIABLE_TYPE(K) == 1 .OR. VARIABLE_TYPE(K) == 2 .OR.  &
                        VARIABLE_TYPE(K) == 3) THEN
                        WRITE (36,"(A,A,A,A,A)") '          ',TRIM(ROUTINE_VARIABLE(K)),  &
                                                 ' = ',TRIM(VARIABLE_INITIALIZATION(K))
                    ENDIF
                ENDIF
             ENDDO

!             If this is the main program, put any multiple precision module variable
!             initializations here.

             IF (ROUTINE_TYPE == 1) THEN
                 REWIND(33)
                 DO
                    READ (33,"(A)",END=110) VAR_LINE
                    IF (VAR_LINE(18:22) == 'Type:') THEN
                        IF (VAR_LINE(24:27) /= '(FM)' .AND. VAR_LINE(24:27) /= '(ZM)' .AND.  &
                            VAR_LINE(24:27) /= '(IM)') CYCLE
                        V_NAME = VAR_LINE(40:71)
                        READ (33,"(A)",END=110) VAR_LINE
                        V_INIT = VAR_LINE(40:132)
                        IF (VAR_LINE(40:132) /= ' ') THEN
                            WRITE (36,"(A,A,A,A,A)") '          ',TRIM(V_NAME),  &
                                                     ' = ',TRIM(V_INIT)
                        ENDIF
                    ENDIF
                 ENDDO
             ENDIF

!             Check for inserting new executable statements at the top of a routine.

  110        IF (JB >= 1) THEN
                 IF (ROUTINE_TYPE /= 1) THEN
                     WRITE (36,"(A)") '      ENDIF'
                 ENDIF
                 WRITE (36,"(A)") ' '
             ELSE
                 IF (ROUTINE_TYPE == 1) THEN
                     WRITE (36,"(A)") ' '
                     WRITE (36,"(A)") '      CALL FM_SET(50)'
                     WRITE (36,"(A)") ' '
                 ELSE IF (ROUTINE_TYPE == 2 .AND. N_VARS > 0) THEN
                     WRITE (36,"(A)") '      CALL FM_ENTER_USER_ROUTINE'
                     WRITE (36,"(A)") ' '
                 ELSE IF (ROUTINE_TYPE == 3) THEN
                     WRITE (36,"(A,A,A)") '      CALL FM_ENTER_USER_FUNCTION('  &
                                          ,TRIM(ROUTINE_NAME),')'
                     WRITE (36,"(A)") ' '
                 ENDIF
             ENDIF
         ENDIF

!             Check for a return statement and add a call to tell the FM routines
!             that this routine is done.

         JB = FIRST_NONBLANK(LINE(J))
         IF (JB > 0) THEN
             TEMP_STATEMENT = LINE(J)
             IF (INDEX_ULC(TEMP_STATEMENT(JB:132),'RETURN') == 1 .AND.  &
                 TEMP_STATEMENT(JB+6:132) == ' ') THEN
                 IF (ROUTINE_TYPE == 2 .AND. N_VARS > 0) THEN
                     WRITE (36,"(A)") '      CALL FM_EXIT_USER_ROUTINE'
                 ELSE IF (ROUTINE_TYPE == 3) THEN
                     WRITE (36,"(A,A,A)") '      CALL FM_EXIT_USER_FUNCTION('  &
                                          ,TRIM(ROUTINE_NAME),')'
                 ENDIF
             ENDIF
         ENDIF

!             Check for an implied return before the end statement of the routine.

         IF (J == LAST .AND. (ROUTINE_TYPE == 2 .OR. ROUTINE_TYPE == 3)) THEN
             JB = LAST_EXECUTABLE(LAST)
             IF (JB > 0) THEN
                 CALL ALL_CAPS(LINE(JB),132,TEMP_STATEMENT)
                 JB = FIRST_NONBLANK(TEMP_STATEMENT)
                 IF (TEMP_STATEMENT(JB:132) /= 'RETURN' .AND. TEMP_STATEMENT(JB:JB+4) /= 'GOTO '  &
                     .AND. TEMP_STATEMENT(JB:JB+5) /= 'GO TO ') THEN
                     IF (ROUTINE_TYPE == 2 .AND. N_VARS > 0) THEN
                         WRITE (36,"(A)") '      CALL FM_EXIT_USER_ROUTINE'
                     ELSE IF (ROUTINE_TYPE == 3) THEN
                         WRITE (36,"(A,A,A)") '      CALL FM_EXIT_USER_FUNCTION('  &
                                          ,TRIM(ROUTINE_NAME),')'
                     ENDIF
                 ENDIF
             ENDIF
         ENDIF

         WRITE (36,"(A)") TRIM(LINE(J))
      ENDDO

      END SUBROUTINE WRITE_PASS2

      SUBROUTINE WRITE_PASS3(LAST)

!  Write the next routine.
!  LINE(1), ..., LINE(LAST) contains the routine.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: END_OF_STATEMENT,START_OF_STATEMENT,J,J1,J2,JB,K,LAST,LAST_NONBLANK,INDEX_ULC

      CALL FIRST_EXECUTABLE(LAST,J1,J2)

      JB = 0
      DO J = 1, LAST
         N_AFTER_READ = 0
         N_AFTER_ALLOC = 0
         IF (J < J2) THEN
             IF (J <= JB) CYCLE
             K = INDEX_ULC(LINE(J),'TYPE (FM)')
             IF (K > 0) THEN
                 IF (LINE(J)(1:MAX(1,K-1)) == ' ' .OR. K == 1) THEN
                     TEMP_STATEMENT = ' '
                     CALL GET_STATEMENT(J,JB,START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT)
                     CALL CONVERT_TYPE_DECLARE(LAST_NONBLANK,END_OF_STATEMENT)
                     CALL WRITE_STATEMENT(LAST_NONBLANK,END_OF_STATEMENT)
                     CYCLE
                 ENDIF
             ENDIF
             K = INDEX_ULC(LINE(J),'TYPE (ZM)')
             IF (K > 0) THEN
                 IF (LINE(J)(1:MAX(1,K-1)) == ' ' .OR. K == 1) THEN
                     TEMP_STATEMENT = ' '
                     CALL GET_STATEMENT(J,JB,START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT)
                     CALL CONVERT_TYPE_DECLARE(LAST_NONBLANK,END_OF_STATEMENT)
                     CALL WRITE_STATEMENT(LAST_NONBLANK,END_OF_STATEMENT)
                     CYCLE
                 ENDIF
             ENDIF
             K = INDEX_ULC(LINE(J),'TYPE (IM)')
             IF (K > 0) THEN
                 IF (LINE(J)(1:MAX(1,K-1)) == ' ' .OR. K == 1) THEN
                     TEMP_STATEMENT = ' '
                     CALL GET_STATEMENT(J,JB,START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT)
                     CALL CONVERT_TYPE_DECLARE(LAST_NONBLANK,END_OF_STATEMENT)
                     CALL WRITE_STATEMENT(LAST_NONBLANK,END_OF_STATEMENT)
                     CYCLE
                 ENDIF
             ENDIF
         ENDIF

!             Add a save statement for any local multiple precision variables.

         IF (J == J2) CALL SAVE_LOCAL_VARS

!             Look for an executable statement where constants should be converted to
!             multiple precision, or where read/write statements are modified.

         IF (J >= J2) THEN
             IF (J <= JB) CYCLE
             CALL GET_STATEMENT(J,JB,START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT)
             TEMP_STATEMENT(LAST_NONBLANK+1:LAST_NONBLANK+1) = ' '
             IF (JB > 0) THEN
                 CALL CONVERT_CONSTANTS(START_OF_STATEMENT,LAST_NONBLANK,END_OF_STATEMENT,1,9999)
                 CALL CONVERT_READS(LAST_NONBLANK,END_OF_STATEMENT)
                 CALL CONVERT_WRITES(LAST_NONBLANK,END_OF_STATEMENT)
             ELSE
                 JB = J
             ENDIF
             CALL WRITE_STATEMENT(LAST_NONBLANK,END_OF_STATEMENT)
             DO K = 1, MAX(N_AFTER_ALLOC,N_AFTER_READ)
                WRITE (23,"(A)") TRIM(AFTER_READ(K))
             ENDDO
             IF (N_AFTER_READ > 0) WRITE (23,"(A)") ' '
             N_AFTER_READ = 0
         ELSE
             WRITE (23,"(A)") TRIM(LINE(J))
             JB = J
         ENDIF

      ENDDO

      END SUBROUTINE WRITE_PASS3

      SUBROUTINE WRITE_STATEMENT(LAST_NONBLANK,END_OF_STATEMENT)

!  TEMP_STATEMENT contains the statement.
!  Write the converted statement to the output file, breaking it into continuation lines
!  if necessary.
!
!  LAST_NONBLANK is the location of the last nonblank character in the statement,
!                including any trailing comments after !.
!  END_OF_STATEMENT is the location of the last nonblank character in the non-comment
!                   part of the statement.

      USE C2FM_VARS
      IMPLICIT NONE
      INTEGER :: END_OF_STATEMENT,JB,JC,K,KA,KB,LAST_NONBLANK
      LOGICAL :: FIRST_LINE

      IF (LAST_NONBLANK < 100) THEN
          WRITE (23,"(A)") TEMP_STATEMENT(1:LAST_NONBLANK)
      ELSE
          KA = 1
          KB = 90
          JC = KB
          FIRST_LINE = .TRUE.
  110     JB = 1
          DO K = KB, KA, -1
             IF (TEMP_STATEMENT(K:K) /= ' ') JB = MAX(K,JB)
             IF (JB > 1 .AND. INDEX(' ,+*/()=<>',TEMP_STATEMENT(K:K)) > 0 .AND.  &
                 QUOTE_LEVEL(K) == 0) THEN
                 JC = K
                 EXIT
             ENDIF
          ENDDO
          IF (JB == LAST_NONBLANK .OR. JB >= END_OF_STATEMENT) THEN
              IF (FIRST_LINE) THEN
                  WRITE (23,"(A)") TRIM(TEMP_STATEMENT(KA:JB))
              ELSE
                  WRITE (23,"(11X,A)") TRIM(TEMP_STATEMENT(KA:JB))
              ENDIF
          ELSE
              IF (FIRST_LINE) THEN
                  WRITE (23,"(A)") TRIM(TEMP_STATEMENT(KA:JC))//'  &'
              ELSE
                  WRITE (23,"(11X,A)") TRIM(TEMP_STATEMENT(KA:JC))//'  &'
              ENDIF
              FIRST_LINE = .FALSE.
              KA = JC + 1
              KB = MIN(KA+78,LAST_NONBLANK)
              GO TO 110
          ENDIF
      ENDIF

      END SUBROUTINE WRITE_STATEMENT
