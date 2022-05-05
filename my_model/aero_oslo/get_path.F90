
subroutine get_aerotab_dir(aerotab_table_dir)

  implicit none

  character (len=255), intent(out)  :: aerotab_table_dir
  character (len=255) :: cwd

  call getcwd(cwd)
  aerotab_table_dir= trim(cwd)//'data/aerotab'

end subroutine get_aerotab_dir
